//! Topology cleanup passes applied after initial polygon extraction.

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

fn polygon_has_unique_boundary_segment<T: GeoFloat>(
    polygons: &[Polygon<T>],
    polygon_index: usize,
) -> bool {
    polygon_unique_boundary_segment_count(polygons, polygon_index) > 0
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

pub(super) fn remove_redundant_overlapping_standalone_polygons<T: GeoFloat>(
    polygons: Vec<Polygon<T>>,
) -> Vec<Polygon<T>> {
    polygons
        .iter()
        .enumerate()
        .filter_map(|(polygon_index, polygon)| {
            if !polygon.interiors().is_empty() {
                return Some(polygon.clone());
            }

            let interior_point = match polygon.interior_point() {
                Some(point) => point,
                None => return Some(polygon.clone()),
            };

            let is_redundant_overlap = polygons.iter().enumerate().any(|(other_index, other)| {
                other_index != polygon_index
                    && other.exterior() != polygon.exterior()
                    && !other.interiors().is_empty()
                    && other.contains(&interior_point)
            });

            if is_redundant_overlap
                && !polygon_has_unique_boundary_segment(&polygons, polygon_index)
            {
                None
            } else {
                Some(polygon.clone())
            }
        })
        .collect()
}

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
