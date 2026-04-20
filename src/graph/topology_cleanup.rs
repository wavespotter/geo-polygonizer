//! Topology cleanup passes applied after initial polygon extraction.
//!
//! The pipeline (defined in [`super::PolygonizerGraph::polygonize`]) runs
//! these passes in order:
//!
//!  1. [`infer_parent_holes_when_output_has_no_holes`] — adopt standalone
//!     polygons as holes of no-hole parents that fully contain them.
//!  2. [`split_touching_boundary_polygons`] — re-polygonize boundaries at
//!     degree>2 nodes to split touching polygons.
//!  3. [`infer_contained_standalone_polygons_as_holes`] — absorb tiny
//!     standalones sitting between existing holes into their parent.
//!  4. [`remove_non_unique_interior_points_for_touching_topology`] — resolve
//!     overlapping polygon ownership by removing the lesser claimant.
//!  5. [`carve_contained_standalones_as_holes`] — carve island polygons as
//!     holes of their enclosing parent (first call).
//!  6. [`merge_touching_holes_in_polygons`] — merge holes connected by
//!     bridge standalones into unified holes.
//!  7. [`carve_contained_standalones_as_holes`] — second carve pass to catch
//!     standalones revealed by merging.

use std::collections::{BTreeMap, BTreeSet};

use crate::nodify::nodify_lines;
use geo::orient::Direction;
use geo::{Area, Contains, GeoFloat, InteriorPoint, Line, LineString, Orient, Polygon, Validation};
use rstar::{Envelope, RTreeObject};

use super::{PolygonizerGraph, ring_to_valid_linestring};

fn polygon_boundary_lines<T: GeoFloat>(polygon: &Polygon<T>) -> Vec<Line<T>> {
    let mut lines: Vec<Line<T>> = polygon.exterior().lines().collect();
    for hole in polygon.interiors() {
        lines.extend(hole.lines());
    }
    lines
}

fn lines_have_same_endpoints<T: GeoFloat>(first: &Line<T>, second: &Line<T>) -> bool {
    (first.start == second.start && first.end == second.end)
        || (first.start == second.end && first.end == second.start)
}

fn polygon_unique_boundary_segment_count<T: GeoFloat>(
    polygons: &[Polygon<T>],
    polygon_index: usize,
) -> usize {
    let polygon_lines = polygon_boundary_lines(&polygons[polygon_index]);

    polygon_lines
        .iter()
        .filter(|line| {
            !polygons
                .iter()
                .enumerate()
                .any(|(other_index, other_polygon)| {
                    if other_index == polygon_index {
                        return false;
                    }

                    polygon_boundary_lines(other_polygon)
                        .iter()
                        .any(|other_line| lines_have_same_endpoints(line, other_line))
                })
        })
        .count()
}

fn polygons_share_any_boundary_segment<T: GeoFloat>(a: &Polygon<T>, b: &Polygon<T>) -> bool {
    let a_lines = polygon_boundary_lines(a);
    let b_lines = polygon_boundary_lines(b);
    a_lines.iter().any(|a_line| {
        b_lines
            .iter()
            .any(|b_line| lines_have_same_endpoints(a_line, b_line))
    })
}

/// When multiple polygons claim the same interior point, remove the lesser
/// claimant.  The "owner" is chosen by most unique boundary segments first,
/// then smallest area.  Polygons that share no boundary with the owner, or
/// that fully contain the owner, are kept.
///
/// Iterates to a fixed point.
pub(super) fn remove_non_unique_interior_points_for_touching_topology<
    T: GeoFloat + rstar::RTreeNum,
>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    let mut current_polygons = polygons;

    loop {
        let mut polygon_indices_to_remove: BTreeSet<usize> = BTreeSet::new();
        let unique_boundary_segment_count_by_index: Vec<usize> = (0..current_polygons.len())
            .map(|polygon_index| {
                polygon_unique_boundary_segment_count(&current_polygons, polygon_index)
            })
            .collect();

        for polygon_index in 0..current_polygons.len() {
            let interior_point = match current_polygons[polygon_index].interior_point() {
                Some(point) => point,
                None => continue,
            };

            let containing_polygon_indices: Vec<usize> = current_polygons
                .iter()
                .enumerate()
                .filter_map(|(candidate_index, candidate_polygon)| {
                    if candidate_polygon.contains(&interior_point) {
                        Some(candidate_index)
                    } else {
                        None
                    }
                })
                .collect();

            if containing_polygon_indices.len() <= 1 {
                continue;
            }

            let owner_index = containing_polygon_indices
                .iter()
                .copied()
                .min_by(|left_index, right_index| {
                    unique_boundary_segment_count_by_index[*right_index]
                        .cmp(&unique_boundary_segment_count_by_index[*left_index])
                        .then(
                            current_polygons[*left_index]
                                .unsigned_area()
                                .total_cmp(&current_polygons[*right_index].unsigned_area()),
                        )
                        .then(left_index.cmp(right_index))
                })
                .expect("owner index should exist when containing polygons are non-empty");

            for candidate_index in containing_polygon_indices {
                if candidate_index != owner_index {
                    let candidate_polygon = &current_polygons[candidate_index];
                    let keep_same_exterior_plain_variant = candidate_polygon.interiors().is_empty()
                        && current_polygons.iter().enumerate().any(
                            |(other_index, other_polygon)| {
                                other_index != candidate_index
                                    && !other_polygon.interiors().is_empty()
                                    && other_polygon.exterior() == candidate_polygon.exterior()
                            },
                        );
                    if keep_same_exterior_plain_variant {
                        continue;
                    }

                    // Don't remove a polygon that doesn't share any boundary
                    // segments with the owner. If they share no boundary, the owner
                    // is simply nested inside the candidate (containment), not a
                    // redundant overlapping polygon.
                    if !polygons_share_any_boundary_segment(
                        candidate_polygon,
                        &current_polygons[owner_index],
                    ) {
                        continue;
                    }

                    // Don't remove a polygon that fully contains the owner.
                    // When the candidate polygon contains the owner, the owner
                    // is a nested face inside the candidate — removing the outer
                    // polygon would destroy a legitimate region.
                    if candidate_polygon.contains(&current_polygons[owner_index]) {
                        continue;
                    }

                    polygon_indices_to_remove.insert(candidate_index);
                }
            }
        }

        if polygon_indices_to_remove.is_empty() {
            break;
        }

        current_polygons = current_polygons
            .into_iter()
            .enumerate()
            .filter_map(|(polygon_index, polygon)| {
                if polygon_indices_to_remove.contains(&polygon_index) {
                    None
                } else {
                    Some(polygon)
                }
            })
            .collect();
    }

    current_polygons
}

/// Re-polygonize each polygon's boundary whenever its boundary graph has
/// degree>2 nodes (i.e. the boundary touches itself).  The polygon is
/// replaced by the interior face sub-polygons.
pub(super) fn split_touching_boundary_polygons<T: GeoFloat + rstar::RTreeNum>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    polygons
        .into_iter()
        .flat_map(split_touching_boundary_polygon)
        .collect()
}

fn build_touching_boundary_graph<T: GeoFloat + rstar::RTreeNum>(
    polygon: &Polygon<T>,
) -> PolygonizerGraph<T> {
    let mut boundary_lines: Vec<_> = polygon.exterior().lines().collect();
    for hole in polygon.interiors() {
        boundary_lines.extend(hole.lines());
    }

    let split_snap_radius =
        T::from(1e-10).unwrap_or(T::epsilon() * T::from(1024.0).unwrap_or(T::one()));
    let noded_boundary_lines = nodify_lines(boundary_lines, split_snap_radius);
    PolygonizerGraph::from_noded_lines(noded_boundary_lines)
}

fn graph_has_touching_topology<T: GeoFloat>(graph: &PolygonizerGraph<T>) -> bool {
    graph
        .nodes_to_outbound_edges
        .values()
        .any(|outbound_edges| outbound_edges.len() > 2)
}

fn extract_face_candidates_inside_polygon<T: GeoFloat>(
    boundary_graph: &PolygonizerGraph<T>,
    container_polygon: &Polygon<T>,
) -> Vec<Polygon<T>> {
    boundary_graph
        .get_minimal_edge_rings()
        .into_iter()
        .filter_map(|ring| {
            let linestring = ring_to_valid_linestring(&ring)?;
            let face_polygon = Polygon::new(linestring, vec![]).orient(Direction::Default);

            let interior_point = face_polygon.interior_point()?;
            if !container_polygon.contains(&interior_point) {
                return None;
            }

            Some(face_polygon)
        })
        .collect()
}

fn deduplicate_faces_by_exterior<T: GeoFloat>(faces: Vec<Polygon<T>>) -> Vec<Polygon<T>> {
    let mut deduplicated_faces: Vec<Polygon<T>> = Vec::new();
    for face in faces {
        let already_present = deduplicated_faces
            .iter()
            .any(|existing| existing.exterior() == face.exterior());
        if !already_present {
            deduplicated_faces.push(face);
        }
    }
    deduplicated_faces
}

fn prune_container_faces<T: GeoFloat>(faces: &[Polygon<T>]) -> Vec<Polygon<T>> {
    faces
        .iter()
        .enumerate()
        .filter_map(|(face_index, face)| {
            let contains_another_face =
                faces
                    .iter()
                    .enumerate()
                    .any(|(other_face_index, other_face)| {
                        if face_index == other_face_index {
                            return false;
                        }

                        let other_interior_point = match other_face.interior_point() {
                            Some(point) => point,
                            None => return false,
                        };

                        face.exterior() != other_face.exterior()
                            && face.contains(&other_interior_point)
                    });

            if contains_another_face {
                None
            } else {
                Some(face.clone())
            }
        })
        .collect()
}

fn select_non_touching_holes<T: GeoFloat + rstar::RTreeNum>(
    polygon: &Polygon<T>,
) -> Vec<LineString<T>> {
    let mut kept_holes: Vec<LineString<T>> = Vec::new();
    for hole in polygon.interiors() {
        let mut candidate_holes = kept_holes.clone();
        candidate_holes.push(hole.clone());

        let candidate_polygon = Polygon::new(polygon.exterior().clone(), candidate_holes.clone())
            .orient(Direction::Default);
        if !graph_has_touching_topology(&build_touching_boundary_graph(&candidate_polygon)) {
            kept_holes = candidate_holes;
        }
    }
    kept_holes
}

fn split_no_hole_polygon_on_repeated_vertex<T: GeoFloat>(
    polygon: &Polygon<T>,
) -> Option<Vec<Polygon<T>>> {
    if !polygon.interiors().is_empty() {
        return None;
    }

    let mut coordinates: Vec<_> = polygon.exterior().points().map(|point| point.0).collect();
    if coordinates.first() == coordinates.last() {
        coordinates.pop();
    }

    if coordinates.len() < 4 {
        return None;
    }

    for first_index in 0..coordinates.len() {
        for second_index in (first_index + 2)..coordinates.len() {
            if first_index == 0 && second_index + 1 == coordinates.len() {
                continue;
            }

            if coordinates[first_index] != coordinates[second_index] {
                continue;
            }

            let first_ring_coords: Vec<_> = coordinates[first_index..=second_index].to_vec();

            let mut second_ring_coords: Vec<_> = coordinates[second_index..].to_vec();
            second_ring_coords.extend_from_slice(&coordinates[..=first_index]);

            if first_ring_coords.len() < 4 || second_ring_coords.len() < 4 {
                continue;
            }

            let first_polygon = Polygon::new(LineString::from(first_ring_coords), vec![])
                .orient(Direction::Default);
            let second_polygon = Polygon::new(LineString::from(second_ring_coords), vec![])
                .orient(Direction::Default);

            if first_polygon.is_valid() && second_polygon.is_valid() {
                return Some(vec![first_polygon, second_polygon]);
            }
        }
    }

    None
}

fn split_touching_boundary_polygon<T: GeoFloat + rstar::RTreeNum>(
    polygon: Polygon<T>,
) -> Vec<Polygon<T>> {
    let boundary_graph = build_touching_boundary_graph(&polygon);
    if !graph_has_touching_topology(&boundary_graph) {
        return vec![polygon];
    }

    let face_candidates = extract_face_candidates_inside_polygon(&boundary_graph, &polygon);
    let deduplicated_faces = deduplicate_faces_by_exterior(face_candidates);
    let split_faces = prune_container_faces(&deduplicated_faces);

    if split_faces.len() >= 2 {
        split_faces
    } else if polygon.interiors().is_empty() {
        split_no_hole_polygon_on_repeated_vertex(&polygon).unwrap_or_else(|| vec![polygon])
    } else {
        let kept_holes = select_non_touching_holes(&polygon);
        vec![Polygon::new(polygon.exterior().clone(), kept_holes).orient(Direction::Default)]
    }
}

/// For each standalone polygon (no holes) that is contained by a parent
/// polygon's exterior, and the parent has no explicitly-assigned holes yet,
/// adopt the standalone as a hole of the smallest qualifying parent.
///
/// This handles the common case where the initial shell/hole assignment
/// did not create a holed polygon because the inner ring was extracted as
/// a standalone shell rather than a hole.
pub(super) fn infer_parent_holes_when_output_has_no_holes<T: GeoFloat + rstar::RTreeNum>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    let polygon_envelopes: Vec<_> = polygons.iter().map(|polygon| polygon.envelope()).collect();

    let mut parent_polygon_index_by_polygon_index: Vec<Option<usize>> = vec![None; polygons.len()];
    for child_polygon_index in 0..polygons.len() {
        let child_envelope = polygon_envelopes[child_polygon_index];
        let mut best_parent: Option<(usize, T)> = None;

        for parent_polygon_index in 0..polygons.len() {
            if child_polygon_index == parent_polygon_index {
                continue;
            }

            let parent_envelope = polygon_envelopes[parent_polygon_index];
            if parent_envelope == child_envelope
                || !parent_envelope.contains_envelope(&child_envelope)
            {
                continue;
            }

            let child_interior_point = match polygons[child_polygon_index].interior_point() {
                Some(point) => point,
                None => continue,
            };

            if !polygons[parent_polygon_index].contains(&child_interior_point) {
                continue;
            }

            let parent_envelope_area = parent_envelope.area();
            if let Some((_, best_area)) = best_parent {
                if parent_envelope_area >= best_area {
                    continue;
                }
            }
            best_parent = Some((parent_polygon_index, parent_envelope_area));
        }

        parent_polygon_index_by_polygon_index[child_polygon_index] =
            best_parent.map(|(parent_polygon_index, _)| parent_polygon_index);
    }

    let mut inferred_holes_by_parent_polygon_index: BTreeMap<usize, Vec<LineString<T>>> =
        BTreeMap::new();
    for (child_polygon_index, parent_polygon_index) in
        parent_polygon_index_by_polygon_index.iter().enumerate()
    {
        let parent_polygon_index = match parent_polygon_index {
            Some(parent_polygon_index) => *parent_polygon_index,
            None => continue,
        };

        let parent_exterior = polygons[parent_polygon_index].exterior();
        let has_same_exterior_polygon_with_explicit_holes =
            polygons.iter().enumerate().any(|(polygon_index, polygon)| {
                polygon_index != parent_polygon_index
                    && polygon.exterior() == parent_exterior
                    && !polygon.interiors().is_empty()
            });
        if has_same_exterior_polygon_with_explicit_holes {
            continue;
        }

        let mut candidate_holes = inferred_holes_by_parent_polygon_index
            .get(&parent_polygon_index)
            .cloned()
            .unwrap_or_default();
        candidate_holes.push(polygons[child_polygon_index].exterior().clone());

        let candidate_polygon = Polygon::new(
            polygons[parent_polygon_index].exterior().clone(),
            candidate_holes.clone(),
        )
        .orient(Direction::Default);
        if !candidate_polygon.is_valid() {
            continue;
        }

        inferred_holes_by_parent_polygon_index.insert(parent_polygon_index, candidate_holes);
    }

    polygons
        .into_iter()
        .enumerate()
        .map(|(polygon_index, polygon)| {
            let mut polygon_holes = polygon.interiors().to_vec();
            polygon_holes.extend(
                inferred_holes_by_parent_polygon_index
                    .get(&polygon_index)
                    .cloned()
                    .unwrap_or_default(),
            );

            Polygon::new(polygon.exterior().clone(), polygon_holes).orient(Direction::Default)
        })
        .collect()
}

/// For each standalone polygon (no holes) whose interior point is contained
/// by another polygon, add the standalone's exterior as a hole of the
/// containing polygon.  This handles "island" polygons that sit fully inside
/// a larger polygon and must be carved out, regardless of whether the parent
/// already has holes.
pub(super) fn carve_contained_standalones_as_holes<T: GeoFloat + rstar::RTreeNum>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    let mut polygons = polygons;

    for child_index in 0..polygons.len() {
        if !polygons[child_index].interiors().is_empty() {
            continue;
        }

        let child_interior_point = match polygons[child_index].interior_point() {
            Some(point) => point,
            None => continue,
        };
        let min_area = T::from(1e-9).unwrap_or(T::epsilon());
        if polygons[child_index].unsigned_area() <= min_area {
            continue;
        }

        let child_exterior = polygons[child_index].exterior().clone();

        // Find the smallest containing polygon that has a different exterior.
        let mut candidate_parents: Vec<(usize, T)> = polygons
            .iter()
            .enumerate()
            .filter_map(|(parent_index, parent_polygon)| {
                if parent_index == child_index {
                    return None;
                }
                if parent_polygon.exterior() == &child_exterior {
                    return None;
                }
                if !parent_polygon.contains(&child_interior_point) {
                    return None;
                }
                Some((parent_index, parent_polygon.unsigned_area()))
            })
            .collect();
        candidate_parents.sort_by(|a, b| a.1.total_cmp(&b.1));

        if candidate_parents.is_empty() {
            continue;
        }

        for (parent_index, _) in candidate_parents {
            // Skip if the child is already a hole of this parent.
            if polygons[parent_index]
                .interiors()
                .iter()
                .any(|hole| hole == &child_exterior)
            {
                continue;
            }

            // Skip if the parent has no holes and a sibling polygon with the
            // same exterior already has holes — the parent is the "plain"
            // variant that should stay hole-free.
            if polygons[parent_index].interiors().is_empty() {
                let sibling_has_holes =
                    polygons
                        .iter()
                        .enumerate()
                        .any(|(other_index, other_polygon)| {
                            other_index != parent_index
                                && other_polygon.exterior() == polygons[parent_index].exterior()
                                && !other_polygon.interiors().is_empty()
                        });
                if sibling_has_holes {
                    continue;
                }
            }

            let mut new_holes = polygons[parent_index].interiors().to_vec();
            new_holes.push(child_exterior.clone());

            let candidate = Polygon::new(polygons[parent_index].exterior().clone(), new_holes)
                .orient(Direction::Default);
            if candidate.is_valid() {
                polygons[parent_index] = candidate;
                break;
            }
        }
    }

    polygons
}

/// Absorb tiny standalone polygons that sit between existing holes into
/// their parent.  Only considers parents with ≥ 2 holes and standalones
/// whose area is below a threshold or that have no unique boundary
/// segments.
pub(super) fn infer_contained_standalone_polygons_as_holes<T: GeoFloat + rstar::RTreeNum>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    let mut polygons = polygons;

    for child_polygon_index in 0..polygons.len() {
        if !polygons[child_polygon_index].interiors().is_empty() {
            continue;
        }

        let child_interior_point = match polygons[child_polygon_index].interior_point() {
            Some(point) => point,
            None => continue,
        };
        let min_child_area = T::from(1e-9).unwrap_or(T::epsilon());
        if polygons[child_polygon_index].unsigned_area() <= min_child_area {
            continue;
        }

        // Only skip standalone polygons with unique boundary segments if
        // their area is above a threshold. Tiny polygons with unique edges
        // are often artifacts (faces between adjacent holes) that should
        // be absorbed as holes of the enclosing polygon.
        let child_area = polygons[child_polygon_index].unsigned_area();
        let tiny_threshold = T::from(0.1).unwrap_or(T::epsilon());
        if child_area > tiny_threshold
            && polygon_unique_boundary_segment_count(&polygons, child_polygon_index) > 0
        {
            continue;
        }

        let child_exterior = polygons[child_polygon_index].exterior().clone();

        let mut candidate_parent_indices: Vec<usize> = polygons
            .iter()
            .enumerate()
            .filter_map(|(parent_polygon_index, parent_polygon)| {
                if parent_polygon_index == child_polygon_index {
                    return None;
                }
                if parent_polygon.interiors().len() < 2 {
                    return None;
                }
                if parent_polygon.exterior() == polygons[child_polygon_index].exterior() {
                    return None;
                }

                // Check containment using only the parent's exterior ring
                // (ignoring existing holes), since the child might be between
                // existing holes.
                let parent_exterior_only = Polygon::new(parent_polygon.exterior().clone(), vec![]);
                if !parent_exterior_only.contains(&child_interior_point) {
                    return None;
                }
                Some(parent_polygon_index)
            })
            .collect();

        candidate_parent_indices.sort_by(|left_index, right_index| {
            polygons[*left_index]
                .unsigned_area()
                .total_cmp(&polygons[*right_index].unsigned_area())
        });

        for parent_polygon_index in candidate_parent_indices {
            let mut candidate_holes = polygons[parent_polygon_index].interiors().to_vec();

            if candidate_holes.iter().any(|hole| hole == &child_exterior) {
                continue;
            }

            candidate_holes.push(child_exterior.clone());
            let candidate_parent_polygon = Polygon::new(
                polygons[parent_polygon_index].exterior().clone(),
                candidate_holes,
            )
            .orient(Direction::Default);
            if !candidate_parent_polygon.is_valid() {
                continue;
            }

            polygons[parent_polygon_index] = candidate_parent_polygon;
            break;
        }
    }

    polygons
}

/// For each polygon with holes, find standalone polygons (from the `polygons`
/// list) that are contained within the polygon's exterior AND share vertices
/// with existing holes. These standalone polygons are "gap fills" between
/// holes, and together with the existing holes they form a connected hole
/// region. Merge all connected hole components (existing holes + matching
/// standalone polygons) into single combined holes.
pub(super) fn merge_touching_holes_in_polygons<T: GeoFloat + rstar::RTreeNum>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    let mut polygons = polygons;

    // For each polygon with >= 2 holes, look for standalone polygons that
    // connect its holes at shared vertices.
    for parent_idx in 0..polygons.len() {
        if polygons[parent_idx].interiors().len() < 2 {
            continue;
        }

        let parent_exterior_only = Polygon::new(polygons[parent_idx].exterior().clone(), vec![]);

        // Build per-hole vertex lists.
        let hole_vertex_sets: Vec<Vec<(T, T)>> = polygons[parent_idx]
            .interiors()
            .iter()
            .map(|hole| hole.0.iter().map(|c| (c.x, c.y)).collect())
            .collect();

        // Collect all hole vertex coords for the parent (flattened).
        let parent_hole_coords: Vec<(T, T)> = hole_vertex_sets
            .iter()
            .flat_map(|s| s.iter().copied())
            .collect();

        // Find standalone polygons (no holes, small) whose exterior shares
        // vertices with the parent's holes and is inside the parent's exterior.
        // Track WHICH holes each standalone connects to.
        let mut connecting_standalone_indices: Vec<usize> = Vec::new();
        let mut any_bridges_multiple_holes = false;
        for child_idx in 0..polygons.len() {
            if child_idx == parent_idx {
                continue;
            }
            if !polygons[child_idx].interiors().is_empty() {
                continue;
            }
            if polygons[child_idx].exterior() == polygons[parent_idx].exterior() {
                continue;
            }
            let child_coords: Vec<(T, T)> = polygons[child_idx]
                .exterior()
                .0
                .iter()
                .map(|c| (c.x, c.y))
                .collect();
            let shares_vertex = child_coords.iter().any(|cv| {
                parent_hole_coords
                    .iter()
                    .any(|pv| pv.0 == cv.0 && pv.1 == cv.1)
            });
            if !shares_vertex {
                continue;
            }
            let child_ip = match polygons[child_idx].interior_point() {
                Some(p) => p,
                None => continue,
            };
            if !parent_exterior_only.contains(&child_ip) {
                continue;
            }
            // Count how many distinct holes this standalone shares vertices with.
            let hole_hits: usize = hole_vertex_sets
                .iter()
                .filter(|hvs| {
                    child_coords
                        .iter()
                        .any(|cv| hvs.iter().any(|hv| hv.0 == cv.0 && hv.1 == cv.1))
                })
                .count();
            if hole_hits >= 2 {
                any_bridges_multiple_holes = true;
            }
            connecting_standalone_indices.push(child_idx);
        }

        // Only merge when at least one standalone bridges 2+ different
        // holes.  If every standalone only touches a single hole, there
        // is nothing to merge and we would lose holes.
        if connecting_standalone_indices.is_empty() || !any_bridges_multiple_holes {
            continue;
        }

        // Collect all edges: parent holes + connecting standalone exteriors.
        let mut all_lines: Vec<Line<T>> = Vec::new();
        for hole in polygons[parent_idx].interiors() {
            all_lines.extend(hole.lines());
        }
        for &child_idx in &connecting_standalone_indices {
            all_lines.extend(polygons[child_idx].exterior().lines());
        }

        // Build graph and extract faces.
        let snap_radius =
            T::from(1e-10).unwrap_or(T::epsilon() * T::from(1024.0).unwrap_or(T::one()));
        let noded_lines = nodify_lines(all_lines, snap_radius);
        let graph = PolygonizerGraph::from_noded_lines(noded_lines);
        let rings = graph.get_minimal_edge_rings();

        // Find the largest ring (outer boundary of the merged holes).
        let mut best_ring: Option<LineString<T>> = None;
        let mut best_area = T::zero();
        for ring in &rings {
            if let Some(ls) = ring_to_valid_linestring(ring) {
                let area = Polygon::new(ls.clone(), vec![]).unsigned_area();
                if area > best_area {
                    best_area = area;
                    best_ring = Some(ls);
                }
            }
        }

        if let Some(merged_ring) = best_ring {
            let candidate =
                Polygon::new(polygons[parent_idx].exterior().clone(), vec![merged_ring])
                    .orient(Direction::Default);
            if candidate.is_valid() {
                polygons[parent_idx] = candidate;
            }
        }
    }

    polygons
}
