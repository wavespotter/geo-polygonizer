use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, btree_map, btree_set};
use std::iter::{FilterMap, Flatten, Map};

use geo::Winding;
use geo::orient::Direction;
use geo::{
    GeoFloat, Line, LineString, LinesIter, MultiPolygon, Orient, Polygon, Validation, coord,
};
use rstar::{Envelope, RTreeObject};

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

#[derive(Debug)]
pub(crate) struct PolygonizerGraph<T: GeoFloat> {
    nodes_to_outbound_edges: BTreeMap<Node<T>, BTreeSet<Edge<T>>>,
}

impl<T: GeoFloat> PolygonizerGraph<T> {
    /// For each directed edge (u -> v), returns the next directed edge along
    /// the face on its left side when traversing the planar graph.
    fn get_edge_to_next_left_face_edge_map(&self) -> BTreeMap<Edge<T>, Edge<T>> {
        let mut next_left_face_edge_by_edge = BTreeMap::new();

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            let outbound_edges: Vec<_> = outbound_edges_at_node.iter().collect();
            let edge_count = outbound_edges.len();
            if edge_count == 0 {
                continue;
            }

            for edge_index in 0..edge_count {
                let reverse_edge = *outbound_edges[edge_index];
                let incoming_edge = reverse_edge.get_symmetrical();
                let next_outgoing_edge =
                    *outbound_edges[(edge_index + edge_count - 1) % edge_count];
                next_left_face_edge_by_edge.insert(incoming_edge, next_outgoing_edge);
            }
        }

        next_left_face_edge_by_edge
    }

    pub(crate) fn from_noded_lines(lines: impl IntoIterator<Item = Line<T>>) -> Self {
        let mut nodes_to_outbound_edges = BTreeMap::new();

        for line in lines.into_iter() {
            let start_node = Node {
                x: line.start.x,
                y: line.start.y,
            };
            let end_node = Node {
                x: line.end.x,
                y: line.end.y,
            };

            nodes_to_outbound_edges
                .entry(start_node)
                .or_insert_with(|| BTreeSet::new())
                .insert(Edge {
                    from: start_node,
                    to: end_node,
                });
            nodes_to_outbound_edges
                .entry(end_node)
                .or_insert_with(|| BTreeSet::new())
                .insert(Edge {
                    from: end_node,
                    to: start_node,
                });
        }

        Self {
            nodes_to_outbound_edges,
        }
    }

    fn get_minimal_edge_rings(&self) -> Vec<Vec<Edge<T>>> {
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

    pub(crate) fn delete_dangles(&mut self) {
        let mut degree_one_nodes: Vec<_> = self
            .nodes_to_outbound_edges
            .iter()
            .filter_map(
                |(node, edges)| {
                    if edges.len() == 1 { Some(*node) } else { None }
                },
            )
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

    pub(crate) fn delete_cut_edges(&mut self) {
        let mut edge_to_index: BTreeMap<Edge<T>, usize> = BTreeMap::new();
        let mut edges_by_index: Vec<Edge<T>> = Vec::new();

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            for edge in outbound_edges_at_node {
                if !edge_to_index.contains_key(edge) {
                    let edge_index = edges_by_index.len();
                    edges_by_index.push(*edge);
                    edge_to_index.insert(*edge, edge_index);
                }
            }
        }

        let mut next_left_face_edge_index_by_edge_index: Vec<Option<usize>> =
            vec![None; edges_by_index.len()];

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            let ordered_outgoing_edges: Vec<_> = outbound_edges_at_node.iter().collect();
            let edge_count = ordered_outgoing_edges.len();
            if edge_count == 0 {
                continue;
            }

            for edge_index in 0..edge_count {
                let reverse_edge = *ordered_outgoing_edges[edge_index];
                let incoming_edge = reverse_edge.get_symmetrical();
                let next_outgoing_edge =
                    *ordered_outgoing_edges[(edge_index + edge_count - 1) % edge_count];

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

        let mut face_label_by_edge_index: Vec<Option<usize>> = vec![None; edges_by_index.len()];
        let mut next_face_label = 0usize;

        for start_edge_index in 0..edges_by_index.len() {
            if face_label_by_edge_index[start_edge_index].is_some() {
                continue;
            }

            let mut traversal_path: Vec<usize> = Vec::new();
            let mut visit_position_by_edge_index: BTreeMap<usize, usize> = BTreeMap::new();
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

    pub(crate) fn polygonize(&self) -> MultiPolygon<T> {
        let edge_rings = self.get_minimal_edge_rings();

        let valid_rings: Vec<_> = edge_rings
            .into_iter()
            .filter_map(|ring| {
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

                // get the linestring back out to return it
                let (linestring, _) = polygon.into_inner();
                Some(linestring)
            })
            .collect();

        let (valid_holes, valid_shells): (Vec<_>, Vec<_>) =
            valid_rings.into_iter().partition(|ring| ring.is_ccw());

        MultiPolygon(
            assign_shells_to_holes(valid_shells, valid_holes)
                .into_iter()
                .filter(|polygon| polygon.is_valid())
                .collect(),
        )
    }
}

type EdgeToLine<'a, T> = fn(&'a Edge<T>) -> Option<Line<T>>;
type EdgeFilterMap<'a, T> = FilterMap<btree_set::Iter<'a, Edge<T>>, EdgeToLine<'a, T>>;
type NodeToEdges<'a, T> = fn((&'a Node<T>, &'a BTreeSet<Edge<T>>)) -> EdgeFilterMap<'a, T>;
type PolygonizerGraphLinesIter<'a, T> =
    Flatten<Map<btree_map::Iter<'a, Node<T>, BTreeSet<Edge<T>>>, NodeToEdges<'a, T>>>;

impl<'a, T: GeoFloat + 'a> LinesIter<'a> for PolygonizerGraph<T> {
    type Scalar = T;
    type Iter = PolygonizerGraphLinesIter<'a, T>;

    fn lines_iter(&'a self) -> Self::Iter {
        self.nodes_to_outbound_edges
            .iter()
            .map(
                (|(_, edges)| {
                    edges.iter().filter_map(
                        (|edge| {
                            if edge.from < edge.to {
                                Some(Line::new(
                                    coord! { x: edge.from.x, y: edge.from.y },
                                    coord! { x: edge.to.x, y: edge.to.y },
                                ))
                            } else {
                                None
                            }
                        }) as EdgeToLine<T>,
                    )
                }) as NodeToEdges<T>,
            )
            .flatten()
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

fn assign_shells_to_holes<T: GeoFloat + rstar::RTreeNum>(
    shells: Vec<LineString<T>>,
    holes: Vec<LineString<T>>,
) -> Vec<Polygon<T>> {
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
        let hole_envelope = hole.envelope();
        let mut matching_shells: Vec<_> = shell_tree
            .locate_in_envelope_intersecting(&hole.envelope())
            .filter(|container| {
                container.envelope.contains_envelope(&hole_envelope)
                    && container.envelope != hole_envelope
            })
            .collect();
        matching_shells.sort_by(|left_shell, right_shell| {
            left_shell
                .envelope
                .area()
                .total_cmp(&right_shell.envelope.area())
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

    for hole_index in assigned_hole_indices {
        let standalone_hole_polygon =
            Polygon::new(holes[hole_index].clone(), vec![]).orient(Direction::Default);
        if !polygons.contains(&standalone_hole_polygon) {
            polygons.push(standalone_hole_polygon);
        }
    }

    polygons
}
