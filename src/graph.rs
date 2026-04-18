use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, btree_map, btree_set};
use std::iter::{FilterMap, Flatten, Map};

use crate::nodify::nodify_lines;
use geo::Winding;
use geo::orient::Direction;
use geo::{
    Area, Contains, GeoFloat, InteriorPoint, Line, LineString, LinesIter, MultiPolygon, Orient,
    Polygon, Validation, coord,
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
        let mut next_left_face_edge_index_by_edge_index: Vec<Option<usize>> = vec![None; edge_count];

        for outbound_edges_at_node in self.nodes_to_outbound_edges.values() {
            let ordered_outgoing_edges: Vec<_> = outbound_edges_at_node.iter().collect();
            let ordered_edge_count = ordered_outgoing_edges.len();
            if ordered_edge_count == 0 {
                continue;
            }

            for edge_index in 0..ordered_edge_count {
                let reverse_edge = *ordered_outgoing_edges[edge_index];
                let incoming_edge = reverse_edge.get_symmetrical();
                let next_outgoing_edge =
                    *ordered_outgoing_edges[(edge_index + ordered_edge_count - 1) % ordered_edge_count];

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
            next_left_face_edge_by_edge.insert(edges_by_index[edge_index], edges_by_index[next_edge_index]);
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
        let (edges_by_index, edge_to_index) = self.get_edges_with_index_map();
        let next_left_face_edge_index_by_edge_index =
            self.get_next_left_face_edge_index_by_edge_index(&edge_to_index, edges_by_index.len());

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

                let (linestring, _) = polygon.into_inner();
                Some(linestring)
            })
            .collect();

        let (valid_holes, valid_shells): (Vec<_>, Vec<_>) =
            valid_rings.into_iter().partition(|ring| ring.is_ccw());

        let valid_polygons: Vec<_> = assign_shells_to_holes(valid_shells, valid_holes)
            .into_iter()
            .filter(|polygon| polygon.is_valid())
            .collect();

        MultiPolygon(remove_redundant_overlapping_standalone_polygons(
            remove_non_unique_interior_points_for_touching_topology(
                split_touching_boundary_polygons(infer_parent_holes_when_output_has_no_holes(
                    valid_polygons,
                )),
            ),
        ))
    }
}

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
            !polygons.iter().enumerate().any(|(other_index, other_polygon)| {
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

fn remove_non_unique_interior_points_for_touching_topology<T: GeoFloat + rstar::RTreeNum>(
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
                        .total_cmp(&current_polygons[*right_index].unsigned_area())
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

fn split_touching_boundary_polygons<T: GeoFloat + rstar::RTreeNum>(
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
            if ring.len() < 3 {
                return None;
            }

            let mut linestring: LineString<T> = ring
                .iter()
                .map(|edge| coord! { x: edge.from.x, y: edge.from.y })
                .collect();
            linestring.close();

            let face_polygon = Polygon::new(linestring, vec![]).orient(Direction::Default);
            if !face_polygon.is_valid() {
                return None;
            }

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
            let contains_another_face = faces.iter().enumerate().any(|(other_face_index, other_face)| {
                if face_index == other_face_index {
                    return false;
                }

                let other_interior_point = match other_face.interior_point() {
                    Some(point) => point,
                    None => return false,
                };

                face.exterior() != other_face.exterior() && face.contains(&other_interior_point)
            });

            if contains_another_face {
                None
            } else {
                Some(face.clone())
            }
        })
        .collect()
}

fn select_non_touching_holes<T: GeoFloat + rstar::RTreeNum>(polygon: &Polygon<T>) -> Vec<LineString<T>> {
    let mut kept_holes: Vec<LineString<T>> = Vec::new();
    for hole in polygon.interiors() {
        let mut candidate_holes = kept_holes.clone();
        candidate_holes.push(hole.clone());

        let candidate_polygon =
            Polygon::new(polygon.exterior().clone(), candidate_holes.clone()).orient(Direction::Default);
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

            let first_polygon = Polygon::new(
                LineString::from(first_ring_coords),
                vec![],
            )
            .orient(Direction::Default);
            let second_polygon = Polygon::new(
                LineString::from(second_ring_coords),
                vec![],
            )
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

fn remove_redundant_overlapping_standalone_polygons<T: GeoFloat>(
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

            if is_redundant_overlap && !polygon_has_unique_boundary_segment(&polygons, polygon_index) {
                None
            } else {
                Some(polygon.clone())
            }
        })
        .collect()
}

fn infer_parent_holes_when_output_has_no_holes<T: GeoFloat + rstar::RTreeNum>(
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
            if parent_envelope == child_envelope || !parent_envelope.contains_envelope(&child_envelope)
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
        let has_same_exterior_polygon_with_explicit_holes = polygons
            .iter()
            .enumerate()
            .any(|(polygon_index, polygon)| {
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

    for hole_index in assigned_hole_indices {
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

    polygons
}


