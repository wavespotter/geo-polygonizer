//! Assigns shell rings to hole rings and builds output polygons.

use std::collections::{BTreeMap, BTreeSet};

use geo::orient::Direction;
use geo::{Area, Contains, GeoFloat, InteriorPoint, LineString, Orient, Polygon, Validation};
use rstar::{Envelope, RTreeObject};

struct ShellContainer<T: GeoFloat> {
    idx: usize,
    envelope: rstar::AABB<geo::Point<T>>,
}

impl<T: GeoFloat + rstar::RTreeNum> RTreeObject for ShellContainer<T> {
    type Envelope = rstar::AABB<geo::Point<T>>;

    fn envelope(&self) -> Self::Envelope {
        self.envelope
    }
}

fn sorted_ring_coords<T: GeoFloat>(ring: &LineString<T>) -> Vec<(T, T)> {
    let mut coords: Vec<(T, T)> = ring.coords().map(|c| (c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    coords.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    coords.dedup();
    coords
}

pub(super) fn assign_shells_to_holes<T: GeoFloat + rstar::RTreeNum>(
    shells: Vec<LineString<T>>,
    holes: Vec<LineString<T>>,
) -> Vec<Polygon<T>> {
    let shell_polygons: Vec<_> = shells
        .iter()
        .cloned()
        .map(|shell| Polygon::new(shell, vec![]))
        .collect();

    let shell_containers = shells
        .iter()
        .enumerate()
        .map(|(idx, shell)| ShellContainer {
            idx,
            envelope: shell.envelope(),
        })
        .collect();
    let shell_tree = rstar::RTree::bulk_load(shell_containers);

    let mut assignments: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

    for (hole_index, hole) in holes.iter().enumerate() {
        let hole_interior_point = match Polygon::new(hole.clone(), vec![]).interior_point() {
            Some(point) => point,
            None => continue,
        };

        let hole_envelope = hole.envelope();
        let mut matching_shells: Vec<_> = shell_tree
            .locate_in_envelope_intersecting(&hole.envelope())
            .filter(|container| {
                container.envelope.contains_envelope(&hole_envelope)
                    && container.envelope != hole_envelope
                    && shell_polygons[container.idx].contains(&hole_interior_point)
            })
            .collect();
        matching_shells.sort_by(|left_shell, right_shell| {
            shell_polygons[left_shell.idx]
                .unsigned_area()
                .total_cmp(&shell_polygons[right_shell.idx].unsigned_area())
        });

        if let Some(container) = matching_shells.first() {
            assignments
                .entry(container.idx)
                .or_insert_with(|| Vec::new())
                .push(hole_index);
        }
    }

    let mut polygons: Vec<Polygon<T>> = shells
        .into_iter()
        .enumerate()
        .map(|(shell_index, shell)| {
            let polygon = Polygon::new(
                shell,
                match assignments.get(&shell_index) {
                    Some(assigned_hole_indices) => assigned_hole_indices
                        .iter()
                        .map(|hole_index| holes[*hole_index].clone())
                        .collect(),
                    None => vec![],
                },
            );
            polygon.orient(Direction::Default)
        })
        .collect();

    let assigned_hole_indices: BTreeSet<usize> = assignments
        .values()
        .flat_map(|hole_indices| hole_indices.iter().copied())
        .collect();

    for hole_index in assigned_hole_indices.iter().copied() {
        let standalone_hole_polygon =
            Polygon::new(holes[hole_index].clone(), vec![]).orient(Direction::Default);

        let overlaps_different_exterior_holed_polygon = polygons.iter().any(|polygon| {
            polygon.exterior() != standalone_hole_polygon.exterior()
                && !polygon.interiors().is_empty()
                && standalone_hole_polygon
                    .exterior()
                    .points()
                    .any(|point| polygon.contains(&point))
        });

        let max_holes_with_same_exterior = polygons
            .iter()
            .filter(|polygon| polygon.exterior() == standalone_hole_polygon.exterior())
            .map(|polygon| polygon.interiors().len())
            .max()
            .unwrap_or(0);

        if !overlaps_different_exterior_holed_polygon
            && max_holes_with_same_exterior <= 1
            && !polygons.contains(&standalone_hole_polygon)
        {
            polygons.push(standalone_hole_polygon);
        }
    }

    // --- Same-envelope hole recovery for invalid shell polygons ---
    //
    // Some hole rings share the exact same bounding box as their
    // containing shell, which causes the `envelope != hole_envelope`
    // filter above to skip them during assignment.  When this happens
    // AND the shell's assigned holes make it invalid (to be filtered
    // later by `polygonize()`), those unassigned same-envelope holes
    // represent genuine face-rings that would otherwise be lost.
    //
    // Recover them as standalone polygons here.  We only do this for
    // shells whose polygon is already invalid, to avoid adding spurious
    // standalone polygons in the common case.

    for (shell_index, assigned) in assignments.iter() {
        if assigned.is_empty() {
            continue;
        }

        // Check whether this shell with its assigned holes is invalid.
        let test_polygon = Polygon::new(
            shell_polygons[*shell_index].exterior().clone(),
            assigned.iter().map(|hi| holes[*hi].clone()).collect(),
        )
        .orient(Direction::Default);

        if test_polygon.is_valid() {
            continue;
        }

        // Shell + holes is invalid — look for unassigned holes that
        // share the shell's bounding box but have different vertices.
        let shell_env = shell_polygons[*shell_index].exterior().envelope();
        let shell_sorted = sorted_ring_coords(shell_polygons[*shell_index].exterior());

        for (hole_index, hole) in holes.iter().enumerate() {
            if assigned_hole_indices.contains(&hole_index) {
                continue;
            }
            if hole.envelope() != shell_env {
                continue;
            }
            let hole_sorted = sorted_ring_coords(hole);
            if hole_sorted == shell_sorted {
                continue; // same vertex set = reverse of shell, skip
            }
            let standalone =
                Polygon::new(hole.clone(), vec![]).orient(Direction::Default);
            if !polygons.contains(&standalone) {
                polygons.push(standalone);
            }
        }
    }

    polygons
}
