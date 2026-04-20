use geo::{GeoFloat, Line, MultiPolygon};

mod graph;
use graph::PolygonizerGraph;

mod nodify;
pub use nodify::nodify_lines;

#[cfg(test)]
mod tests;

/// Polygonizes an input line set into a [`geo::MultiPolygon`].
///
/// Input lines are expected to already be suitably noded if intersections or
/// overlaps should be split at intersection points. Use [`nodify_lines`] first
/// when this precondition is not guaranteed.
pub fn polygonize<T: GeoFloat>(lines: impl IntoIterator<Item = Line<T>>) -> MultiPolygon<T> {
    let mut graph = PolygonizerGraph::from_noded_lines(lines);
    graph.delete_dangles();
    graph.delete_cut_edges();

    graph.polygonize()
}
