//! LinesIter implementation for the polygonizer graph.

use std::collections::{BTreeSet, btree_map, btree_set};
use std::iter::{FilterMap, Flatten, Map};

use geo::{GeoFloat, Line, LinesIter, coord};

use super::{Edge, Node, PolygonizerGraph};

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
