//! Planar graph representation and polygon extraction pipeline.
//!
//! This module stores noded linework as directed edges sorted in CCW order per
//! origin node, then traverses left faces to assemble minimal rings.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use geo::Winding;
use geo::{GeoFloat, Line, LineString, MultiPolygon, Polygon, Validation, coord};

mod lines_iter;
mod shell_assignment;
mod topology_cleanup;

use shell_assignment::assign_shells_to_holes;

type EdgeIndex = usize;
type FaceLabel = usize;

#[derive(Copy, Clone, Debug)]
pub(crate) struct Node<T: GeoFloat> {
    x: T,
    y: T,
}

impl<T: GeoFloat> PartialEq for Node<T> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}
impl<T: GeoFloat> Eq for Node<T> {}
impl<T: GeoFloat> PartialOrd for Node<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<T: GeoFloat> Ord for Node<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.x.total_cmp(&other.x).then(self.y.total_cmp(&other.y))
    }
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct Edge<T: GeoFloat> {
    from: Node<T>,
    to: Node<T>,
}

impl<T: GeoFloat> PartialEq for Edge<T> {
    fn eq(&self, other: &Self) -> bool {
        self.from == other.from && self.to == other.to
    }
}
impl<T: GeoFloat> Eq for Edge<T> {}
impl<T: GeoFloat> PartialOrd for Edge<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
// sort edges in CCW order, assuming same from-node
impl<T: GeoFloat> Ord for Edge<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.from.cmp(&other.from) {
            Ordering::Equal => {
                let self_dx = self.to.x - self.from.x;
                let self_dy = self.to.y - self.from.y;
                let other_dx = other.to.x - other.from.x;
                let other_dy = other.to.y - other.from.y;

                let self_is_upper_half =
                    self_dy > T::zero() || (self_dy == T::zero() && self_dx >= T::zero());
                let other_is_upper_half =
                    other_dy > T::zero() || (other_dy == T::zero() && other_dx >= T::zero());

                match self_is_upper_half.cmp(&other_is_upper_half).reverse() {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => {
                        let cross_product = self_dx * other_dy - self_dy * other_dx;
                        if cross_product > T::zero() {
                            return Ordering::Less;
                        }
                        if cross_product < T::zero() {
                            return Ordering::Greater;
                        }

                        let self_length_squared = self_dx * self_dx + self_dy * self_dy;
                        let other_length_squared = other_dx * other_dx + other_dy * other_dy;
                        self_length_squared
                            .total_cmp(&other_length_squared)
                            .then(self.to.cmp(&other.to))
                    }
                }
            }
            ordering => ordering,
        }
    }
}

impl<T: GeoFloat> Edge<T> {
    fn get_symmetrical(&self) -> Self {
        Edge {
            from: self.to,
            to: self.from,
        }
    }
}

fn line_endpoints_to_nodes<T: GeoFloat>(line: Line<T>) -> (Node<T>, Node<T>) {
    (
        Node {
            x: line.start.x,
            y: line.start.y,
        },
        Node {
            x: line.end.x,
            y: line.end.y,
        },
    )
}

fn ring_to_valid_linestring<T: GeoFloat>(ring: &[Edge<T>]) -> Option<LineString<T>> {
    if ring.len() < 3 {
        return None;
    }

    let mut linestring: LineString<T> = ring
        .iter()
        .map(|edge| coord! { x: edge.from.x, y: edge.from.y })
        .collect();
    linestring.close();
    let polygon = Polygon::new(linestring, Default::default());
    if !polygon.is_valid() {
        return None;
    }

    let (linestring, _) = polygon.into_inner();
    Some(linestring)
}

/// Splits an edge ring at any repeated "from" node (pinch point).
///
/// When the left-face traversal visits the same node twice in one ring
/// (a figure-8 or multi-lobed shape connected by single-point pinch
/// vertices), this function breaks the ring into sub-rings, each of
/// which forms a simple loop through the pinch node exactly once.
/// The result is applied recursively so that rings with multiple pinch
/// points are fully decomposed.
fn split_edge_ring_at_pinch_points<T: GeoFloat>(ring: Vec<Edge<T>>) -> Vec<Vec<Edge<T>>> {
    // Build a map: from-node → first position where we saw it.
    let mut first_occurrence: BTreeMap<Node<T>, usize> = BTreeMap::new();
    let mut split_at: Option<(usize, usize)> = None;

    for (position, edge) in ring.iter().enumerate() {
        if let Some(&first_position) = first_occurrence.get(&edge.from) {
            split_at = Some((first_position, position));
            break;
        }
        first_occurrence.insert(edge.from, position);
    }

    let (first_position, second_position) = match split_at {
        Some(positions) => positions,
        None => return vec![ring], // no repeated node, nothing to split
    };

    // ring[first_position..second_position] forms one sub-ring (the loop
    // that starts and ends at the pinch node).
    let inner_ring: Vec<Edge<T>> = ring[first_position..second_position].to_vec();

    // The remainder is everything before + everything from second_position onward.
    let mut outer_ring: Vec<Edge<T>> = ring[..first_position].to_vec();
    outer_ring.extend_from_slice(&ring[second_position..]);

    // Recurse on both halves in case there are additional pinch points.
    let mut result = Vec::new();
    result.extend(split_edge_ring_at_pinch_points(inner_ring));
    result.extend(split_edge_ring_at_pinch_points(outer_ring));
    result
}

#[derive(Debug)]
pub(crate) struct PolygonizerGraph<T: GeoFloat> {
    nodes_to_outbound_edges: BTreeMap<Node<T>, BTreeSet<Edge<T>>>,
}

impl<T: GeoFloat> PolygonizerGraph<T> {
    fn insert_undirected_edge(
        nodes_to_outbound_edges: &mut BTreeMap<Node<T>, BTreeSet<Edge<T>>>,
        start_node: Node<T>,
        end_node: Node<T>,
    ) {
        nodes_to_outbound_edges
            .entry(start_node)
            .or_insert_with(BTreeSet::new)
            .insert(Edge {
                from: start_node,
                to: end_node,
            });
        nodes_to_outbound_edges
            .entry(end_node)
            .or_insert_with(BTreeSet::new)
            .insert(Edge {
                from: end_node,
                to: start_node,
            });
    }

    fn compute_face_label_by_edge_index(
        edge_count: usize,
        next_left_face_edge_index_by_edge_index: &[Option<EdgeIndex>],
    ) -> Vec<Option<FaceLabel>> {
        let mut face_label_by_edge_index: Vec<Option<FaceLabel>> = vec![None; edge_count];
        let mut next_face_label: FaceLabel = 0;

        for start_edge_index in 0..edge_count {
            if face_label_by_edge_index[start_edge_index].is_some() {
                continue;
            }

            let mut traversal_path: Vec<EdgeIndex> = Vec::new();
            let mut visit_position_by_edge_index: BTreeMap<EdgeIndex, usize> = BTreeMap::new();
            let mut current_edge_index = start_edge_index;

            loop {
                if let Some(existing_face_label) = face_label_by_edge_index[current_edge_index] {
                    for traversed_edge_index in traversal_path {
                        face_label_by_edge_index[traversed_edge_index] = Some(existing_face_label);
                    }
                    break;
                }

                if let Some(&cycle_start_position) =
                    visit_position_by_edge_index.get(&current_edge_index)
                {
                    let face_label = next_face_label;
                    next_face_label += 1;

                    for &traversed_edge_index in &traversal_path[cycle_start_position..] {
                        face_label_by_edge_index[traversed_edge_index] = Some(face_label);
                    }
                    for &traversed_edge_index in &traversal_path[..cycle_start_position] {
                        face_label_by_edge_index[traversed_edge_index] = Some(face_label);
                    }
                    break;
                }

                visit_position_by_edge_index.insert(current_edge_index, traversal_path.len());
                traversal_path.push(current_edge_index);

                current_edge_index =
                    match next_left_face_edge_index_by_edge_index[current_edge_index] {
                        Some(next_edge_index) => next_edge_index,
                        None => break,
                    };
            }
        }

        face_label_by_edge_index
    }

    fn collect_undirected_cut_edges(
        edges_by_index: &[Edge<T>],
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        face_label_by_edge_index: &[Option<FaceLabel>],
    ) -> BTreeSet<(Node<T>, Node<T>)> {
        let mut undirected_edges_to_remove: BTreeSet<(Node<T>, Node<T>)> = BTreeSet::new();
        for (edge_index, edge) in edges_by_index.iter().enumerate() {
            if edge.from > edge.to {
                continue;
            }

            let symmetric_edge_index = *edge_to_index
                .get(&edge.get_symmetrical())
                .expect("symmetric edge index should exist");

            if face_label_by_edge_index[edge_index]
                == face_label_by_edge_index[symmetric_edge_index]
            {
                undirected_edges_to_remove.insert((edge.from, edge.to));
            }
        }
        undirected_edges_to_remove
    }

    fn remove_undirected_edges(
        &mut self,
        undirected_edges_to_remove: &BTreeSet<(Node<T>, Node<T>)>,
    ) {
        self.nodes_to_outbound_edges
            .retain(|_, outbound_edges_at_node| {
                outbound_edges_at_node.retain(|edge| {
                    let undirected_edge = if edge.from <= edge.to {
                        (edge.from, edge.to)
                    } else {
                        (edge.to, edge.from)
                    };
                    !undirected_edges_to_remove.contains(&undirected_edge)
                });
                !outbound_edges_at_node.is_empty()
            });
    }

    fn get_edges_with_index_map(&self) -> (Vec<Edge<T>>, BTreeMap<Edge<T>, usize>) {
        let mut edge_to_index: BTreeMap<Edge<T>, usize> = BTreeMap::new();
        let mut edges_by_index: Vec<Edge<T>> = Vec::new();

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            for edge in outbound_edges_at_node {
                if edge_to_index.contains_key(edge) {
                    continue;
                }

                let edge_index = edges_by_index.len();
                edges_by_index.push(*edge);
                edge_to_index.insert(*edge, edge_index);
            }
        }

        (edges_by_index, edge_to_index)
    }

    fn get_next_left_face_edge_index_by_edge_index(
        &self,
        edge_to_index: &BTreeMap<Edge<T>, usize>,
        edge_count: usize,
    ) -> Vec<Option<usize>> {
        let mut next_left_face_edge_index_by_edge_index: Vec<Option<usize>> =
            vec![None; edge_count];

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            let ordered_outgoing_edges: Vec<_> = outbound_edges_at_node.iter().collect();
            let ordered_edge_count = ordered_outgoing_edges.len();
            if ordered_edge_count == 0 {
                continue;
            }

            for edge_index in 0..ordered_edge_count {
                let reverse_edge = *ordered_outgoing_edges[edge_index];
                let incoming_edge = reverse_edge.get_symmetrical();
                let next_outgoing_edge = *ordered_outgoing_edges
                    [(edge_index + ordered_edge_count - 1) % ordered_edge_count];

                let incoming_edge_index = *edge_to_index
                    .get(&incoming_edge)
                    .expect("incoming edge index should exist");
                let next_outgoing_edge_index = *edge_to_index
                    .get(&next_outgoing_edge)
                    .expect("next outgoing edge index should exist");

                next_left_face_edge_index_by_edge_index[incoming_edge_index] =
                    Some(next_outgoing_edge_index);
            }
        }

        next_left_face_edge_index_by_edge_index
    }

    /// For each directed edge (u -> v), returns the next directed edge along
    /// the face on its left side when traversing the planar graph.
    fn get_edge_to_next_left_face_edge_map(&self) -> BTreeMap<Edge<T>, Edge<T>> {
        let (edges_by_index, edge_to_index) = self.get_edges_with_index_map();
        let next_left_face_edge_index_by_edge_index =
            self.get_next_left_face_edge_index_by_edge_index(&edge_to_index, edges_by_index.len());

        let mut next_left_face_edge_by_edge = BTreeMap::new();
        for (edge_index, next_edge_index) in next_left_face_edge_index_by_edge_index
            .iter()
            .copied()
            .enumerate()
        {
            let Some(next_edge_index) = next_edge_index else {
                continue;
            };
            next_left_face_edge_by_edge
                .insert(edges_by_index[edge_index], edges_by_index[next_edge_index]);
        }

        next_left_face_edge_by_edge
    }

    /// Builds a directed graph from already-noded linework.
    ///
    /// Each undirected segment is represented as two directed edges with opposite directions.
    pub(crate) fn from_noded_lines(lines: impl IntoIterator<Item = Line<T>>) -> Self {
        let mut nodes_to_outbound_edges = BTreeMap::new();

        for line in lines.into_iter() {
            let (start_node, end_node) = line_endpoints_to_nodes(line);
            Self::insert_undirected_edge(&mut nodes_to_outbound_edges, start_node, end_node);
        }

        Self {
            nodes_to_outbound_edges,
        }
    }

    pub(super) fn get_minimal_edge_rings(&self) -> Vec<Vec<Edge<T>>> {
        let next_left_face_edge_by_edge = self.get_edge_to_next_left_face_edge_map();

        let mut rings = Vec::new();
        let mut visited_directed_edges = BTreeSet::new();
        for edge in next_left_face_edge_by_edge.keys() {
            if visited_directed_edges.contains(edge) {
                continue;
            }

            let mut ring = vec![*edge];
            visited_directed_edges.insert(*edge);
            let mut next_edge = *edge;
            loop {
                next_edge = match next_left_face_edge_by_edge.get(&next_edge) {
                    Some(next) => *next,
                    _ => break,
                };
                if next_edge == *edge {
                    break;
                }
                if visited_directed_edges.contains(&next_edge) {
                    break;
                }
                ring.push(next_edge);
                visited_directed_edges.insert(next_edge);
            }
            rings.push(ring);
        }
        rings
    }

    /// Iteratively removes dangling degree-1 chains from the graph.
    pub(crate) fn delete_dangles(&mut self) {
        let mut degree_one_nodes: Vec<_> = self
            .nodes_to_outbound_edges
            .iter()
            .filter_map(|(node, edges)| (edges.len() == 1).then_some(*node))
            .collect();

        loop {
            let source_node = match degree_one_nodes.pop() {
                Some(node) => node,
                _ => break,
            };

            let mut source_edge_list = match self.nodes_to_outbound_edges.remove(&source_node) {
                Some(source_edge_list) => source_edge_list,
                _ => continue,
            };

            let source_edge = match source_edge_list.pop_last() {
                Some(source_edge) => source_edge,
                _ => continue,
            };

            // remove the symmetric edge
            let dest_node = source_edge.to;
            let dest_edge_list = match self.nodes_to_outbound_edges.get_mut(&dest_node) {
                Some(dest_edge_list) => dest_edge_list,
                _ => continue,
            };
            dest_edge_list.remove(&source_edge.get_symmetrical());

            // if the target node is now degree-one (or zero), add to the stack
            if dest_edge_list.len() <= 1 {
                degree_one_nodes.push(dest_node);
            }
        }
    }

    /// Removes cut edges by labeling directed-edge faces and dropping edges whose
    /// two directions belong to the same face.
    pub(crate) fn delete_cut_edges(&mut self) {
        let (edges_by_index, edge_to_index) = self.get_edges_with_index_map();
        let next_left_face_edge_index_by_edge_index =
            self.get_next_left_face_edge_index_by_edge_index(&edge_to_index, edges_by_index.len());

        let face_label_by_edge_index = Self::compute_face_label_by_edge_index(
            edges_by_index.len(),
            &next_left_face_edge_index_by_edge_index,
        );
        let undirected_edges_to_remove = Self::collect_undirected_cut_edges(
            &edges_by_index,
            &edge_to_index,
            &face_label_by_edge_index,
        );
        self.remove_undirected_edges(&undirected_edges_to_remove);
    }

    /// Extracts valid rings, classifies shell/hole orientation, and applies
    /// downstream topology cleanup passes.
    pub(crate) fn polygonize(&self) -> MultiPolygon<T>
    where
        T: rstar::RTreeNum,
    {
        let edge_rings = self.get_minimal_edge_rings();

        // Split any figure-8 rings at pinch-point nodes before the validity
        // filter so that both sub-regions survive as separate polygons.
        let split_rings: Vec<_> = edge_rings
            .into_iter()
            .flat_map(split_edge_ring_at_pinch_points)
            .collect();

        let valid_rings: Vec<_> = split_rings
            .into_iter()
            .filter_map(|ring| ring_to_valid_linestring(&ring))
            .collect();

        let (valid_holes, valid_shells): (Vec<_>, Vec<_>) =
            valid_rings.into_iter().partition(|ring| ring.is_ccw());

        let valid_polygons: Vec<_> = assign_shells_to_holes(valid_shells, valid_holes)
            .into_iter()
            .filter(|polygon| polygon.is_valid())
            .collect();

        // --- Topology cleanup pipeline (7 passes) ---
        // Assign standalone polygons as holes of no-hole parents that contain them.
        let polygons =
            topology_cleanup::infer_parent_holes_when_output_has_no_holes(valid_polygons);
        // Re-polygonize boundaries at degree>2 nodes to split touching polygons.
        let polygons =
            topology_cleanup::split_touching_boundary_polygons(polygons);
        // Absorb tiny standalones sitting between holes into their parent polygon.
        let polygons =
            topology_cleanup::infer_contained_standalone_polygons_as_holes(polygons);
        // Resolve overlapping polygon ownership by removing shared interior points.
        let polygons =
            topology_cleanup::remove_non_unique_interior_points_for_touching_topology(polygons);
        // Carve contained standalone polygons as holes of their enclosing parent.
        let polygons =
            topology_cleanup::carve_contained_standalones_as_holes(polygons);
        // Merge holes connected by bridge standalones into unified holes.
        let polygons =
            topology_cleanup::merge_touching_holes_in_polygons(polygons);
        // Second carve pass to catch standalones revealed by merging.
        let polygons =
            topology_cleanup::carve_contained_standalones_as_holes(polygons);

        MultiPolygon(polygons)
    }
}
