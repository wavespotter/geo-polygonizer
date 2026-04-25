//! Planar graph representation and polygon extraction pipeline.
//!
//! This module stores noded linework as directed edges sorted in CCW order per
//! origin node, then follows the JTS polygonizer ring traversal:
//! `computeNextCWEdges`, `findLabeledEdgeRings`, and
//! `convertMaximalToMinimalEdgeRings`.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use geo::Winding;
use geo::{GeoFloat, Line, LineString, MultiPolygon, Polygon, Validation, coord};

mod lines_iter;
mod shell_assignment;

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

    fn collect_undirected_cut_edges(
        edges_by_index: &[Edge<T>],
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        edge_ring_label_by_edge_index: &[Option<FaceLabel>],
    ) -> BTreeSet<(Node<T>, Node<T>)> {
        let mut undirected_edges_to_remove: BTreeSet<(Node<T>, Node<T>)> = BTreeSet::new();
        for (edge_index, edge) in edges_by_index.iter().enumerate() {
            if edge.from > edge.to {
                continue;
            }

            let symmetric_edge_index = *edge_to_index
                .get(&edge.get_symmetrical())
                .expect("symmetric edge index should exist");

            if edge_ring_label_by_edge_index[edge_index].is_some()
                && edge_ring_label_by_edge_index[edge_index]
                    == edge_ring_label_by_edge_index[symmetric_edge_index]
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

    /// Analogous to JTS `computeNextCWEdges`.
    ///
    /// For every directed edge entering a node, choose the next outgoing edge
    /// in clockwise order around that node.
    fn compute_next_cw_edges_by_edge_index(
        &self,
        edge_to_index: &BTreeMap<Edge<T>, usize>,
        edge_count: usize,
    ) -> Vec<Option<usize>> {
        let mut next_cw_edge_index_by_edge_index: Vec<Option<usize>> = vec![None; edge_count];

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            let ordered_outgoing_edges: Vec<_> = outbound_edges_at_node.iter().copied().collect();
            let ordered_edge_count = ordered_outgoing_edges.len();
            if ordered_edge_count == 0 {
                continue;
            }

            for edge_index in 0..ordered_edge_count {
                let outgoing_edge = ordered_outgoing_edges[edge_index];
                let incoming_edge = outgoing_edge.get_symmetrical();
                let next_outgoing_edge = ordered_outgoing_edges
                    [(edge_index + ordered_edge_count - 1) % ordered_edge_count];

                let incoming_edge_index = *edge_to_index
                    .get(&incoming_edge)
                    .expect("incoming edge index should exist");
                let next_outgoing_edge_index = *edge_to_index
                    .get(&next_outgoing_edge)
                    .expect("next outgoing edge index should exist");

                next_cw_edge_index_by_edge_index[incoming_edge_index] =
                    Some(next_outgoing_edge_index);
            }
        }

        next_cw_edge_index_by_edge_index
    }

    /// Analogous to JTS `findDirEdgesInRing`.
    fn find_dir_edges_in_ring(
        start_edge_index: EdgeIndex,
        next_edge_index_by_edge_index: &[Option<EdgeIndex>],
    ) -> Vec<EdgeIndex> {
        let mut ring = Vec::new();
        let mut visited_edge_indices = BTreeSet::new();
        let mut current_edge_index = start_edge_index;

        loop {
            if !visited_edge_indices.insert(current_edge_index) {
                break;
            }

            ring.push(current_edge_index);

            let Some(next_edge_index) = next_edge_index_by_edge_index[current_edge_index] else {
                break;
            };
            if next_edge_index == start_edge_index {
                break;
            }

            current_edge_index = next_edge_index;
        }

        ring
    }

    /// Analogous to JTS `findLabeledEdgeRings`.
    fn find_labeled_edge_rings(
        edge_count: usize,
        next_edge_index_by_edge_index: &[Option<EdgeIndex>],
    ) -> (Vec<Option<FaceLabel>>, Vec<EdgeIndex>) {
        let mut edge_ring_label_by_edge_index = vec![None; edge_count];
        let mut maximal_ring_start_edge_indices = Vec::new();
        let mut next_edge_ring_label: FaceLabel = 0;

        for start_edge_index in 0..edge_count {
            if edge_ring_label_by_edge_index[start_edge_index].is_some() {
                continue;
            }

            maximal_ring_start_edge_indices.push(start_edge_index);
            let ring_edge_indices =
                Self::find_dir_edges_in_ring(start_edge_index, next_edge_index_by_edge_index);
            for edge_index in ring_edge_indices {
                edge_ring_label_by_edge_index[edge_index] = Some(next_edge_ring_label);
            }
            next_edge_ring_label += 1;
        }

        (
            edge_ring_label_by_edge_index,
            maximal_ring_start_edge_indices,
        )
    }

    fn get_degree_for_label(
        &self,
        node: Node<T>,
        label: FaceLabel,
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        edge_ring_label_by_edge_index: &[Option<FaceLabel>],
    ) -> usize {
        self.nodes_to_outbound_edges
            .get(&node)
            .map(|outbound_edges_at_node| {
                outbound_edges_at_node
                    .iter()
                    .filter(|edge| {
                        let edge_index = *edge_to_index
                            .get(edge)
                            .expect("edge index should exist for graph edge");
                        edge_ring_label_by_edge_index[edge_index] == Some(label)
                    })
                    .count()
            })
            .unwrap_or(0)
    }

    /// Analogous to JTS `findIntersectionNodes`.
    fn find_intersection_nodes(
        &self,
        start_edge_index: EdgeIndex,
        label: FaceLabel,
        edges_by_index: &[Edge<T>],
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        edge_ring_label_by_edge_index: &[Option<FaceLabel>],
        next_edge_index_by_edge_index: &[Option<EdgeIndex>],
    ) -> BTreeSet<Node<T>> {
        let mut intersection_nodes = BTreeSet::new();

        for edge_index in
            Self::find_dir_edges_in_ring(start_edge_index, next_edge_index_by_edge_index)
        {
            let node = edges_by_index[edge_index].from;
            if self.get_degree_for_label(node, label, edge_to_index, edge_ring_label_by_edge_index)
                > 1
            {
                intersection_nodes.insert(node);
            }
        }

        intersection_nodes
    }

    /// Analogous to JTS `computeNextCCWEdges`.
    ///
    /// This rewires only the directed edges that have the supplied maximal-ring
    /// label, which is what subdivides maximal edge rings into minimal rings.
    fn compute_next_ccw_edges_at_node(
        &self,
        node: Node<T>,
        label: FaceLabel,
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        edge_ring_label_by_edge_index: &[Option<FaceLabel>],
        next_edge_index_by_edge_index: &mut [Option<EdgeIndex>],
    ) {
        let Some(outbound_edges_at_node) = self.nodes_to_outbound_edges.get(&node) else {
            return;
        };

        let ordered_outgoing_edges: Vec<_> = outbound_edges_at_node.iter().copied().collect();
        let mut first_out_edge_index = None;
        let mut previous_in_edge_index = None;

        for outgoing_edge in ordered_outgoing_edges.iter().rev() {
            let outgoing_edge_index = *edge_to_index
                .get(outgoing_edge)
                .expect("outgoing edge index should exist");
            let incoming_edge_index = *edge_to_index
                .get(&outgoing_edge.get_symmetrical())
                .expect("incoming edge index should exist");

            let out_edge_index = (edge_ring_label_by_edge_index[outgoing_edge_index]
                == Some(label))
            .then_some(outgoing_edge_index);
            let in_edge_index = (edge_ring_label_by_edge_index[incoming_edge_index] == Some(label))
                .then_some(incoming_edge_index);

            if out_edge_index.is_none() && in_edge_index.is_none() {
                continue;
            }

            if let Some(in_edge_index) = in_edge_index {
                previous_in_edge_index = Some(in_edge_index);
            }

            if let Some(out_edge_index) = out_edge_index {
                if let Some(in_edge_index) = previous_in_edge_index.take() {
                    next_edge_index_by_edge_index[in_edge_index] = Some(out_edge_index);
                }
                if first_out_edge_index.is_none() {
                    first_out_edge_index = Some(out_edge_index);
                }
            }
        }

        if let Some(in_edge_index) = previous_in_edge_index {
            if let Some(out_edge_index) = first_out_edge_index {
                next_edge_index_by_edge_index[in_edge_index] = Some(out_edge_index);
            }
        }
    }

    /// Analogous to JTS `convertMaximalToMinimalEdgeRings`.
    fn convert_maximal_to_minimal_edge_rings(
        &self,
        maximal_ring_start_edge_indices: &[EdgeIndex],
        edges_by_index: &[Edge<T>],
        edge_to_index: &BTreeMap<Edge<T>, EdgeIndex>,
        edge_ring_label_by_edge_index: &[Option<FaceLabel>],
        next_edge_index_by_edge_index: &mut [Option<EdgeIndex>],
    ) {
        for &start_edge_index in maximal_ring_start_edge_indices {
            let Some(label) = edge_ring_label_by_edge_index[start_edge_index] else {
                continue;
            };

            let intersection_nodes = self.find_intersection_nodes(
                start_edge_index,
                label,
                edges_by_index,
                edge_to_index,
                edge_ring_label_by_edge_index,
                next_edge_index_by_edge_index,
            );

            for node in intersection_nodes {
                self.compute_next_ccw_edges_at_node(
                    node,
                    label,
                    edge_to_index,
                    edge_ring_label_by_edge_index,
                    next_edge_index_by_edge_index,
                );
            }
        }
    }

    /// Analogous to JTS `getEdgeRings`.
    pub(super) fn get_edge_rings(&self) -> Vec<Vec<Edge<T>>> {
        let (edges_by_index, edge_to_index) = self.get_edges_with_index_map();
        let mut next_edge_index_by_edge_index =
            self.compute_next_cw_edges_by_edge_index(&edge_to_index, edges_by_index.len());

        let (edge_ring_label_by_edge_index, maximal_ring_start_edge_indices) =
            Self::find_labeled_edge_rings(edges_by_index.len(), &next_edge_index_by_edge_index);

        self.convert_maximal_to_minimal_edge_rings(
            &maximal_ring_start_edge_indices,
            &edges_by_index,
            &edge_to_index,
            &edge_ring_label_by_edge_index,
            &mut next_edge_index_by_edge_index,
        );

        let mut rings = Vec::new();
        let mut visited_edge_indices = BTreeSet::new();
        for start_edge_index in 0..edges_by_index.len() {
            if visited_edge_indices.contains(&start_edge_index) {
                continue;
            }

            let mut ring = Vec::new();
            let mut current_edge_index = start_edge_index;

            loop {
                if !visited_edge_indices.insert(current_edge_index) {
                    break;
                }

                ring.push(edges_by_index[current_edge_index]);

                let Some(next_edge_index) = next_edge_index_by_edge_index[current_edge_index]
                else {
                    break;
                };
                if next_edge_index == start_edge_index {
                    break;
                }

                current_edge_index = next_edge_index;
            }

            if !ring.is_empty() {
                rings.push(ring);
            }
        }

        rings
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
        let next_cw_edge_index_by_edge_index =
            self.compute_next_cw_edges_by_edge_index(&edge_to_index, edges_by_index.len());

        let (edge_ring_label_by_edge_index, _) =
            Self::find_labeled_edge_rings(edges_by_index.len(), &next_cw_edge_index_by_edge_index);
        let undirected_edges_to_remove = Self::collect_undirected_cut_edges(
            &edges_by_index,
            &edge_to_index,
            &edge_ring_label_by_edge_index,
        );
        self.remove_undirected_edges(&undirected_edges_to_remove);
    }

    /// Extracts valid rings, classifies shell/hole orientation, assigns
    /// holes to shells (analogous to JTS `polygonize`), and applies
    /// optional topology cleanup for self-touching boundaries.
    ///
    /// Shell/hole classification follows geo's OGC convention:
    ///   CCW ring = shell (bounded face exterior)
    ///   CW  ring = hole  (reverse-face / interior ring)
    ///
    /// This matches JTS's `computeHole()` adapted for OGC winding:
    /// JTS uses CW-exterior convention; geo-rust uses CCW-exterior,
    /// so CCW=shell here.
    pub(crate) fn polygonize(&self) -> MultiPolygon<T>
    where
        T: rstar::RTreeNum,
    {
        let edge_rings = self.get_edge_rings();

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

        // Analogous to JTS `findShellsAndHoles` — classify by winding.
        // In geo's OGC convention: CCW = shell (exterior), CW = hole.
        let (valid_shells, valid_holes): (Vec<_>, Vec<_>) =
            valid_rings.into_iter().partition(|ring| ring.is_ccw());

        // Analogous to JTS `assignHolesToShells` + `extractPolygons`.
        let valid_polygons: Vec<_> = assign_shells_to_holes(valid_shells, valid_holes)
            .into_iter()
            .filter(|polygon| polygon.is_valid())
            .collect();

        MultiPolygon(valid_polygons)
    }
}
