//! Like JTS `HoleAssigner` — assigns hole rings to the smallest containing shell.
//!
//! Each hole (CW ring in OGC convention) is assigned to the smallest shell
//! (CCW ring) whose polygon strictly contains the hole's interior point.
//! Unassigned holes (e.g. the CW reverse-face of each shell's own boundary)
//! are silently dropped, matching JTS behaviour.

use std::collections::BTreeMap;

use geo::Winding;
use geo::orient::Direction;
use geo::{Contains, GeoFloat, InteriorPoint, LineString, Orient, Polygon};
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

/// Like JTS `assignHolesToShells` + `extractPolygons`.
///
/// For each hole, finds the smallest containing shell (by area) whose
/// envelope strictly contains the hole's envelope and whose polygon
/// contains the hole's interior point. Builds one `Polygon` per shell
/// (shell exterior + assigned holes), oriented per OGC convention.
pub(super) fn assign_shells_to_holes<T: GeoFloat + rstar::RTreeNum>(
    shells: Vec<LineString<T>>,
    holes: Vec<LineString<T>>,
) -> Vec<Polygon<T>> {
    let shell_polygons: Vec<_> = shells
        .iter()
        .cloned()
        .map(|shell| Polygon::new(shell, vec![]))
        .collect();

    let shell_containers: Vec<_> = shells
        .iter()
        .enumerate()
        .map(|(idx, shell)| ShellContainer {
            idx,
            envelope: shell.envelope(),
        })
        .collect();
    let shell_tree = rstar::RTree::bulk_load(shell_containers);

    let mut assignments: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

    // Analogous to JTS `assignHoleToShell` — for each hole, find the
    // smallest containing shell (analogous to JTS `findEdgeRingContaining`).
    for (hole_index, hole) in holes.iter().enumerate() {
        // Hole rings are CW. To get a reliable interior point inside the
        // geometric shape, orient the ring CCW before computing.
        let mut oriented_hole = hole.clone();
        if !oriented_hole.is_ccw() {
            oriented_hole.make_ccw_winding();
        }
        let hole_interior_point = match Polygon::new(oriented_hole, vec![]).interior_point() {
            Some(point) => point,
            None => continue,
        };

        let hole_envelope = hole.envelope();

        // Envelope must strictly contain the hole envelope (not equal).
        let mut matching_shells: Vec<_> = shell_tree
            .locate_in_envelope_intersecting(&hole_envelope)
            .filter(|container| {
                container.envelope.contains_envelope(&hole_envelope)
                    && container.envelope != hole_envelope
                    && shell_polygons[container.idx].contains(&hole_interior_point)
            })
            .collect();

        // Pick smallest containing shell by envelope area (O(1) per shell,
        // vs O(n) for polygon area). Equivalent for the simple closed rings
        // the polygonizer produces, matching JTS's approach.
        matching_shells.sort_by(|a, b| a.envelope.area().total_cmp(&b.envelope.area()));

        if let Some(container) = matching_shells.first() {
            assignments
                .entry(container.idx)
                .or_default()
                .push(hole_index);
        }
        // Unassigned holes are silently dropped (JTS behaviour).
    }

    // Analogous to JTS `extractPolygons` — build one Polygon per shell.
    let result: Vec<_> = shells
        .into_iter()
        .enumerate()
        .map(|(shell_index, shell)| {
            let assigned_holes: Vec<LineString<T>> = assignments
                .get(&shell_index)
                .map(|indices| indices.iter().map(|&i| holes[i].clone()).collect())
                .unwrap_or_default();

            Polygon::new(shell, assigned_holes).orient(Direction::Default)
        })
        .collect();
    result
}
