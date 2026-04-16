use geo::{GeoFloat, Line, MultiPolygon};

mod graph;
use graph::PolygonizerGraph;

mod nodify;
pub use nodify::nodify_lines;

#[cfg(test)]
mod tests;

pub fn polygonize<T: GeoFloat>(lines: impl IntoIterator<Item = Line<T>>) -> MultiPolygon<T> {
    let mut graph = PolygonizerGraph::from_noded_lines(lines);
    graph.delete_dangles();
    graph.delete_cut_edges();

    graph.polygonize()
}
