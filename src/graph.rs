use std::cmp::Ordering;
use std::iter::{FilterMap, Flatten, Map};
use std::collections::{BTreeMap, BTreeSet, btree_map, btree_set};

use geo::orient::Direction;
use geo::{GeoFloat, Line, LineString, LinesIter, MultiPolygon, Orient, Polygon, Validation, Winding, coord};
use rstar::{Envelope, RTreeObject};

#[derive(PartialOrd, Ord, PartialEq, Eq)]
enum Quadrant {
    NE = 0,
    NW = 1,
    SW = 2,
    SE = 3,
}

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
                match self.quadrant().cmp(&other.quadrant()) {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => {
                        let dx1 = self.to.x - self.from.x;
                        let dy1 = self.to.y - self.from.y;
                        let dx2 = other.to.x - other.from.x;
                        let dy2 = other.to.y - other.from.y;

                        let cross_cmp = (dx1 * dy2 - dx2 * dy1).total_cmp(&T::zero());
                        if cross_cmp != Ordering::Equal {
                            return cross_cmp;
                        }

                        let len1 = dx1 * dx1 + dy1 * dy1;
                        let len2 = dx2 * dx2 + dy2 * dy2;
                        len1.total_cmp(&len2).then(self.to.cmp(&other.to))
                    }
                }
            }
            ordering => ordering,
        }
    }
}

impl<T: GeoFloat> Edge<T> {
    fn quadrant(&self) -> Quadrant {
        let dx = self.to.x - self.from.x;
        let dy = self.to.y - self.from.y;

        match (dx >= T::zero(), dy >= T::zero()) {
            (true, true) => Quadrant::NE,
            (true, false) => Quadrant::SE,
            (false, true) => Quadrant::NW,
            (false, false) => Quadrant::SW,
        }
    }

    fn get_symmetrical(&self) -> Self {
        Edge { from: self.to, to: self.from }
    }
}

#[derive(Debug)]
pub(crate) struct PolygonizerGraph<T: GeoFloat> {
    nodes_to_outbound_edges: BTreeMap<Node<T>, BTreeSet<Edge<T>>>,
}

impl<T: GeoFloat> PolygonizerGraph<T> {
    pub(crate) fn from_noded_lines(lines: impl IntoIterator<Item = Line<T>>) -> Self {
        let mut nodes_to_outbound_edges = BTreeMap::new();

        for line in lines.into_iter() {
            let node1 = Node { x: line.start.x, y: line.start.y };
            let node2 = Node { x: line.end.x, y: line.end.y };

            nodes_to_outbound_edges
                .entry(node1)
                .or_insert_with(|| BTreeSet::new())
                .insert(Edge { from: node1, to: node2 });
            nodes_to_outbound_edges
                .entry(node2)
                .or_insert_with(|| BTreeSet::new())
                .insert(Edge { from: node2, to: node1 });
        }

        Self { nodes_to_outbound_edges }
    }

    /// For each directed edge (from → to) in the graph, returns a map to the
    /// next outbound edge from `from` in clockwise order (wrapping cyclically).
    fn get_edge_to_next_cw_edge_map(&self) -> BTreeMap<Edge<T>, Edge<T>> {
        let mut map = BTreeMap::new();
        for (_node, edges) in &self.nodes_to_outbound_edges {
            let edges_vec: Vec<_> = edges.iter().collect();
            let n = edges_vec.len();
            for i in 0..n {
                let current_edge = edges_vec[i];
                let prev_edge = edges_vec[(i + n - 1) % n];
                map.insert(*current_edge, prev_edge.get_symmetrical());
            }
        }
        map
    }

    /// Labels all edges by which maximal ring they belong to, and returns the
    /// start edge for each ring (one entry per distinct ring label).
    fn label_edges_by_ring(&self, edge_map: &BTreeMap<Edge<T>, Edge<T>>) -> (BTreeMap<Edge<T>, usize>, Vec<Edge<T>>) {
        let mut labels: BTreeMap<Edge<T>, usize> = BTreeMap::new();
        let mut ring_starts: Vec<Edge<T>> = Vec::new();
        let mut label: usize = 0;
        for edge in edge_map.keys() {
            if labels.contains_key(edge) {
                continue;
            }

            ring_starts.push(*edge);
            labels.insert(*edge, label);
            let mut next = *edge;
            loop {
                next = match edge_map.get(&next) {
                    Some(next) => *next,
                    _ => break
                };
                if next == *edge {
                    break;
                }
                if labels.contains_key(&next) {
                    break;
                }
                labels.insert(next, label);
            }
            label += 1;
        }
        (labels, ring_starts)
    }

    /// Converts the maximal edge ring map (CW traversal) into a minimal edge ring map.
    ///
    /// Mirrors JTS `convertMaximalToMinimalEdgeRings`: for each maximal ring, traverse
    /// it to find "intersection nodes" (nodes where that ring has degree > 1), then
    /// rewire only those nodes for that ring's label. This avoids incorrectly rewiring
    /// nodes that are only visited once per ring.
    fn convert_maximal_to_minimal_edge_rings(
        &self,
        edge_map: &BTreeMap<Edge<T>, Edge<T>>,
        ring_labeling: &BTreeMap<Edge<T>, usize>,
        ring_starts: &[Edge<T>],
    ) -> BTreeMap<Edge<T>, Edge<T>> {
        let mut map = edge_map.clone();

        for &start_edge in ring_starts {
            let label = match ring_labeling.get(&start_edge) {
                Some(&l) => l,
                None => continue,
            };

            // Traverse this maximal ring via CW pointers to find intersection nodes:
            // nodes where the ring's label has degree > 1 (i.e. the ring visits the
            // node more than once, meaning the maximal ring self-intersects there).
            let mut intersection_nodes: Vec<Node<T>> = Vec::new();
            let mut cur = start_edge;
            loop {
                let node = cur.from;
                // Count how many edges at this node belong to this ring label
                // (both outgoing and incoming directions).
                // Count only outgoing edges at this node that carry this ring's label.
                // This mirrors JTS `getDegree(node, label)`, which counts only outgoing
                // directed edges. A node has degree > 1 only when the ring visits it more
                // than once (i.e., the maximal ring self-intersects there). Counting
                // incoming edges too would make every normal traversal node look like an
                // intersection (outgoing 1 + incoming 1 = 2) and cause all nodes to be
                // spuriously rewired.
                let degree = self.nodes_to_outbound_edges
                    .get(&node)
                    .map(|edges| edges.iter().filter(|e| {
                        ring_labeling.get(*e) == Some(&label)
                    }).count())
                    .unwrap_or(0);
                if degree > 1 {
                    intersection_nodes.push(node);
                }
                cur = match edge_map.get(&cur) {
                    Some(&next) => next,
                    None => break,
                };
                if cur == start_edge {
                    break;
                }
            }

            // Rewire the next-pointers at each intersection node for this ring label,
            // turning the maximal ring into minimal rings (CCW traversal order).
            for node in intersection_nodes {
                let edges = match self.nodes_to_outbound_edges.get(&node) {
                    Some(e) => e,
                    None => continue,
                };
                let edges_vec: Vec<&Edge<T>> = edges.iter().collect();
                let n = edges_vec.len();

                let mut first_out_edge: Option<Edge<T>> = None;
                let mut prev_in_edge: Option<Edge<T>> = None;

                // Iterate in reverse (CW) order — same as JTS iterating CW around the star.
                for i in (0..n).rev() {
                    let edge = edges_vec[i];
                    let symmetrical = edge.get_symmetrical();

                    let out_edge = if ring_labeling.get(edge) == Some(&label) { Some(*edge) } else { None };
                    let in_edge = if ring_labeling.get(&symmetrical) == Some(&label) { Some(symmetrical) } else { None };

                    if out_edge.is_none() && in_edge.is_none() {
                        continue;
                    }

                    if let Some(in_edge) = in_edge {
                        prev_in_edge = Some(in_edge);
                    }

                    if let Some(out_edge) = out_edge {
                        if let Some(prev_in) = prev_in_edge.take() {
                            map.insert(prev_in, out_edge);
                        }
                        if first_out_edge.is_none() {
                            first_out_edge = Some(out_edge);
                        }
                    }
                }

                // Wrap-around: connect the last unmatched in_edge to the first out_edge.
                if let (Some(prev_in), Some(first_out)) = (prev_in_edge, first_out_edge) {
                    map.insert(prev_in, first_out);
                }
            }
        }

        map
    }

    fn get_minimal_edge_rings(&self) -> Vec<Vec<Edge<T>>> {
        let cw_map = self.get_edge_to_next_cw_edge_map();
        let (ring_labels, ring_starts) = self.label_edges_by_ring(&cw_map);
        let minimal_map = self.convert_maximal_to_minimal_edge_rings(&cw_map, &ring_labels, &ring_starts);

        let mut rings = Vec::new();
        let mut is_in_ring = BTreeSet::new();
        for edge in minimal_map.keys() {
            if is_in_ring.contains(edge) {
                continue;
            }

            let mut ring = vec![*edge];
            is_in_ring.insert(*edge);
            let mut next = *edge;
            loop {
                next = match minimal_map.get(&next) {
                    Some(next) => *next,
                    _ => break
                };
                if next == *edge {
                    break;
                }
                if is_in_ring.contains(&next) {
                    break;
                }
                ring.push(next);
                is_in_ring.insert(next);
            }
            rings.push(ring);
        }
        rings
    }

    pub(crate) fn delete_dangles(&mut self) {
        let mut degree_one_nodes: Vec<_> = self.nodes_to_outbound_edges.iter().filter_map(|(node, edges)| {
            if edges.len() == 1 {
                Some(*node)
            } else {
                None
            }
        }).collect();

        loop {
            let source_node = match degree_one_nodes.pop() {
                Some(node) => node,
                _ => break
            };

            let mut source_edge_list = match self.nodes_to_outbound_edges.remove(&source_node) {
                Some(source_edge_list) => source_edge_list,
                _ => continue
            };

            let source_edge = match source_edge_list.pop_last() {
                Some(source_edge) => source_edge,
                _ => continue
            };

            // remove the symmetric edge
            let dest_node = source_edge.to;
            let dest_edge_list = match self.nodes_to_outbound_edges.get_mut(&dest_node) {
                Some(dest_edge_list) => dest_edge_list,
                _ => continue
            };
            dest_edge_list.remove(&source_edge.get_symmetrical());

            // if the target node is now degree-one (or zero), add to the stack
            if dest_edge_list.len() <= 1 {
                degree_one_nodes.push(dest_node);
            }
        }
    }

    pub(crate) fn delete_cut_edges(&mut self) {
        let cw_map = self.get_edge_to_next_cw_edge_map();
        let (ring_labels, _) = self.label_edges_by_ring(&cw_map);
        for (edge, label) in ring_labels.iter() {
            if ring_labels.get(&edge.get_symmetrical()) == Some(label) {
                self.nodes_to_outbound_edges.get_mut(&edge.from).map(|edges| edges.remove(edge));
                self.nodes_to_outbound_edges.get_mut(&edge.to).map(|edges| edges.remove(&edge.get_symmetrical()));
            }
        }
    }

    pub(crate) fn polygonize(&self) -> MultiPolygon<T> {
        let edge_rings = self.get_minimal_edge_rings();

        let (valid_holes, valid_shells): (Vec<_>, Vec<_>) = edge_rings.into_iter().filter_map(|ring| {
            if ring.len() < 3 {
                return None;
            }

            let mut linestring: LineString<T> = ring.iter().map(|edge| coord! { x: edge.from.x, y: edge.from.y }).collect();
            linestring.close();
            let polygon = Polygon::new(linestring, Default::default());
            if !polygon.is_valid() {
                return None;
            }

            // get the linestring back out to return it
            let (linestring, _) = polygon.into_inner();
            Some(linestring)
        }).partition(|linestring| linestring.is_ccw());

        MultiPolygon(assign_shells_to_holes(valid_shells, valid_holes))
    }
}

type EdgeToLine<'a, T> = fn(&'a Edge<T>) -> Option<Line<T>>;
type EdgeFilterMap<'a, T> = FilterMap<btree_set::Iter<'a, Edge<T>>, EdgeToLine<'a, T>>;
type NodeToEdges<'a, T> = fn((&'a Node<T>, &'a BTreeSet<Edge<T>>)) -> EdgeFilterMap<'a, T>;
type PolygonizerGraphLinesIter<'a, T> = Flatten<Map<btree_map::Iter<'a, Node<T>, BTreeSet<Edge<T>>>, NodeToEdges<'a, T>>>;

impl<'a, T: GeoFloat + 'a> LinesIter<'a> for PolygonizerGraph<T> {
    type Scalar = T;
    type Iter = PolygonizerGraphLinesIter<'a, T>;

    fn lines_iter(&'a self) -> Self::Iter {
        self.nodes_to_outbound_edges.iter().map(
            (|(_, edges)| {
                edges.iter().filter_map(
                    (|edge| {
                        if edge.from < edge.to {
                            Some(Line::new(coord! { x: edge.from.x, y: edge.from.y }, coord! { x: edge.to.x, y: edge.to.y }))
                        } else {
                            None
                        }
                    }) as EdgeToLine<T>
                )
            }) as NodeToEdges<T>
        ).flatten()
    }
}

struct ShellContainer<T: GeoFloat> {
    idx: usize,
    envelope: rstar::AABB<geo::Point<T>>,
}

impl<T: GeoFloat + rstar::RTreeNum> rstar::RTreeObject for ShellContainer<T> {
    type Envelope = rstar::AABB<geo::Point<T>>;

    fn envelope(&self) -> Self::Envelope {
        self.envelope
    }
}

fn assign_shells_to_holes<T: GeoFloat + rstar::RTreeNum>(shells: Vec<LineString<T>>, holes: Vec<LineString<T>>) -> Vec<Polygon<T>> {
    let shell_containers = shells.iter().enumerate().map(|(idx, shell)| ShellContainer { idx, envelope: shell.envelope() }).collect();
    let shell_tree = rstar::RTree::bulk_load(shell_containers);

    let mut assignments: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

    for (i, hole) in holes.iter().enumerate() {
        let hole_envelope = hole.envelope();
        let mut matching: Vec<_> = shell_tree.locate_in_envelope_intersecting(&hole.envelope()).filter(|container| container.envelope.contains_envelope(&hole_envelope) && container.envelope != hole_envelope).collect();
        matching.sort_by(|c1, c2| c1.envelope.area().total_cmp(&c2.envelope.area()));

        if let Some(container) = matching.first() {
            assignments
                .entry(container.idx)
                .or_insert_with(|| Vec::new())
                .push(i);
        }
    }

    shells.into_iter().enumerate().map(|(i, shell)| {
        let polygon = Polygon::new(
            shell,
            match assignments.get(&i) {
                Some(assignments) => assignments.iter().map(|j| holes[*j].clone()).collect(),
                None => vec![]
            }
        );
        polygon.orient(Direction::Default)
    }).collect()
}
