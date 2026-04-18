use std::fs;
use std::path::PathBuf;

use geo::{Line, LinesIter};
use geo_polygonizer::polygonize;
use geojson::{Feature, FeatureCollection, Geometry, GeometryValue};

fn fixture_path(kind: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("polygonizer")
        .join(kind)
        .join(format!("{name}.geojson"))
}

fn load_input_lines(name: &str) -> Vec<Line<f64>> {
    let path = fixture_path("input", name);
    let text = fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read fixture {}: {err}", path.display()));
    let collection: FeatureCollection = serde_json::from_str(&text)
        .unwrap_or_else(|err| panic!("failed to parse fixture {}: {err}", path.display()));

    collection
        .features
        .into_iter()
        .flat_map(|feature| {
            let geometry = feature
                .geometry
                .unwrap_or_else(|| panic!("input fixture {name} has a feature with no geometry"));
            match geometry.value {
                GeometryValue::LineString { .. } => {
                    let linestring: geo::LineString<f64> = geometry.try_into().unwrap();
                    linestring.lines().collect::<Vec<_>>()
                }
                GeometryValue::MultiLineString { .. } => {
                    let multiline: geo::MultiLineString<f64> = geometry.try_into().unwrap();
                    multiline.lines_iter().collect::<Vec<_>>()
                }
                other => {
                    panic!("input fixture {name} uses unsupported geometry type: {other:?}")
                }
            }
        })
        .collect()
}

fn main() {
    let name = std::env::args().nth(1).unwrap_or_else(|| "minimal_secondary_probe_split".to_string());
    let input_lines = load_input_lines(&name);
    let polygons = polygonize(input_lines);

    let features = polygons
        .0
        .into_iter()
        .map(|polygon| {
            let geometry = Geometry::new(GeometryValue::from(&polygon));
            Feature {
                bbox: None,
                geometry: Some(geometry),
                id: None,
                properties: None,
                foreign_members: None,
            }
        })
        .collect();

    let collection = FeatureCollection {
        bbox: None,
        features,
        foreign_members: None,
    };

    let output_path = fixture_path("output", &name);
    let output_json = serde_json::to_string_pretty(&collection).expect("serialize output");
    fs::write(&output_path, output_json)
        .unwrap_or_else(|err| panic!("failed to write fixture {}: {err}", output_path.display()));
}
