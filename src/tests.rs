use geo::{LinesIter, MultiLineString};

use super::*;

#[test]
fn nodify_preserves_first_seen_coordinates_when_not_collapsed() {
    let input = vec![geo::Line::new(
        geo::Coord { x: 0.14_f64, y: 0.23_f64 },
        geo::Coord { x: 1.27_f64, y: 2.31_f64 },
    )];

    let output: Vec<_> = nodify::nodify_lines(input, 0.1_f64);

    assert_eq!(output.len(), 1);
    let line = output[0];

    let forward = line.start.x == 0.14_f64
        && line.start.y == 0.23_f64
        && line.end.x == 1.27_f64
        && line.end.y == 2.31_f64;
    let backward = line.start.x == 1.27_f64
        && line.start.y == 2.31_f64
        && line.end.x == 0.14_f64
        && line.end.y == 0.23_f64;

    assert!(
        forward || backward,
        "expected exact cached endpoints, got {:?}",
        line
    );
}

#[test]
fn nodify_collapsed_points_use_first_seen_coordinate() {
    let first = geo::Coord { x: 0.14_f64, y: 0.0_f64 };
    let second = geo::Coord { x: 0.149_f64, y: 0.0_f64 };

    let input = vec![
        geo::Line::new(first, geo::Coord { x: 0.5_f64, y: 0.0_f64 }),
        geo::Line::new(second, geo::Coord { x: 0.9_f64, y: 0.0_f64 }),
    ];

    let output: Vec<_> = nodify::nodify_lines(input, 0.1_f64);

    let mut saw_collapsed_key = false;
    for line in output {
        for p in [line.start, line.end] {
            let qx = (p.x / 0.1_f64).round() as i64;
            let qy = (p.y / 0.1_f64).round() as i64;
            if qx == 1 && qy == 0 {
                saw_collapsed_key = true;
                assert_eq!(p.x, first.x, "collapsed x should use first-seen value");
                assert_eq!(p.y, first.y, "collapsed y should use first-seen value");
            }
        }
    }

    assert!(saw_collapsed_key, "expected at least one endpoint on the collapsed key");
}

#[test]
fn lines_work() {
    let data: geojson::FeatureCollection = serde_json::from_str(include_str!("complex.geojson")).unwrap();
    let data: Vec<geo::LineString<f64>> = data.features.iter().map(|feat| feat.geometry.as_ref().unwrap().try_into().unwrap()).collect();
    let data = geo::MultiLineString(data);

    let mut graph = PolygonizerGraph::from_noded_lines(data.lines_iter());
    graph.delete_dangles();
    graph.delete_cut_edges();

    let polygons = graph.polygonize();

    let polygons: geojson::Geometry = (&polygons).try_into().unwrap();
    let polygons = polygons.to_string();
    println!("{polygons}");
}

#[test]
fn holes_work() {
    let data: geojson::FeatureCollection = serde_json::from_str(include_str!("holes.geojson")).unwrap();
    let data: Vec<geo::Polygon<f64>> = data.features.iter().map(|feat| feat.geometry.as_ref().unwrap().try_into().unwrap()).collect();
    let data = geo::MultiPolygon(data);

    let mut graph = PolygonizerGraph::from_noded_lines(data.lines_iter());
    graph.delete_dangles();
    graph.delete_cut_edges();

    let polygons = graph.polygonize();

    let polygons: geojson::Geometry = (&polygons).try_into().unwrap();
    let polygons = polygons.to_string();
    println!("{polygons}");
}

#[test]
fn nodify_and_polygonize() {
    let data: geojson::FeatureCollection = serde_json::from_str(include_str!("crossing.geojson")).unwrap();
    let data: Vec<geo::LineString<f64>> = data.features.iter().map(|feat| feat.geometry.as_ref().unwrap().try_into().unwrap()).collect();
    let data = geo::MultiLineString(data);

    let nodified = nodify_lines(data.lines_iter(), 1e-9);

    let mls: MultiLineString = nodified.iter().cloned().collect();

    let mut graph = PolygonizerGraph::from_noded_lines(nodified.into_iter());
    graph.delete_dangles();
    graph.delete_cut_edges();

    let polygons = graph.polygonize();

    let polygons: geojson::Geometry = (&polygons).try_into().unwrap();
    let polygons = polygons.to_string();
    println!("{polygons}");

    let polygons: geojson::Geometry = (&mls).try_into().unwrap();
    let polygons = polygons.to_string();
    println!("{polygons}");
}