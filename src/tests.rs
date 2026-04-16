use std::fs;
use std::path::PathBuf;

use geo::{Line, LinesIter, MultiPolygon, Polygon};
use geojson::GeometryValue;

use super::*;

#[test]
fn nodify_preserves_first_seen_coordinates_when_not_collapsed() {
    let input = vec![geo::Line::new(
        geo::Coord {
            x: 0.14_f64,
            y: 0.23_f64,
        },
        geo::Coord {
            x: 1.27_f64,
            y: 2.31_f64,
        },
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
    let first = geo::Coord {
        x: 0.14_f64,
        y: 0.0_f64,
    };
    let second = geo::Coord {
        x: 0.149_f64,
        y: 0.0_f64,
    };

    let input = vec![
        geo::Line::new(
            first,
            geo::Coord {
                x: 0.5_f64,
                y: 0.0_f64,
            },
        ),
        geo::Line::new(
            second,
            geo::Coord {
                x: 0.9_f64,
                y: 0.0_f64,
            },
        ),
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

    assert!(
        saw_collapsed_key,
        "expected at least one endpoint on the collapsed key"
    );
}

fn fixture_path(kind: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("polygonizer")
        .join(kind)
        .join(format!("{name}.geojson"))
}

fn load_feature_collection(kind: &str, name: &str) -> geojson::FeatureCollection {
    let path = fixture_path(kind, name);
    let text = fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read fixture {}: {err}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|err| panic!("failed to parse fixture {}: {err}", path.display()))
}

fn load_input_lines(name: &str) -> Vec<Line<f64>> {
    let collection = load_feature_collection("input", name);
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

fn load_expected_polygons_from(output_kind: &str, name: &str) -> MultiPolygon<f64> {
    let collection = load_feature_collection(output_kind, name);
    let mut polygons = Vec::new();

    for feature in collection.features {
        let geometry = feature
            .geometry
            .unwrap_or_else(|| panic!("output fixture {name} has a feature with no geometry"));
        match geometry.value {
            GeometryValue::Polygon { .. } => {
                let polygon: Polygon<f64> = geometry.try_into().unwrap();
                polygons.push(polygon);
            }
            GeometryValue::MultiPolygon { .. } => {
                let multipolygon: MultiPolygon<f64> = geometry.try_into().unwrap();
                polygons.extend(multipolygon.0);
            }
            other => {
                panic!("output fixture {name} uses unsupported geometry type: {other:?}");
            }
        }
    }

    MultiPolygon(polygons)
}

fn canonicalize_ring(ring: &geo::LineString<f64>) -> Vec<(i64, i64)> {
    let mut coords: Vec<(i64, i64)> = ring
        .points()
        .map(|point| (point.x().round() as i64, point.y().round() as i64))
        .collect();

    if coords.first() == coords.last() {
        coords.pop();
    }

    fn remove_collinear_vertices(coords: &[(i64, i64)]) -> Vec<(i64, i64)> {
        if coords.len() <= 3 {
            return coords.to_vec();
        }

        let len = coords.len();
        let mut out = Vec::with_capacity(len);

        for i in 0..len {
            let prev = coords[(i + len - 1) % len];
            let cur = coords[i];
            let next = coords[(i + 1) % len];

            let previous_to_current_dx = cur.0 - prev.0;
            let previous_to_current_dy = cur.1 - prev.1;
            let current_to_next_dx = next.0 - cur.0;
            let current_to_next_dy = next.1 - cur.1;
            let cross_product = previous_to_current_dx * current_to_next_dy
                - previous_to_current_dy * current_to_next_dx;

            if cross_product != 0 {
                out.push(cur);
            }
        }

        if out.len() < 3 { coords.to_vec() } else { out }
    }

    coords = remove_collinear_vertices(&coords);

    fn best_rotation(mut ring: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
        let len = ring.len();
        if len <= 1 {
            return ring;
        }

        let mut best = ring.clone();
        for _ in 1..len {
            ring.rotate_left(1);
            if ring < best {
                best = ring.clone();
            }
        }
        best
    }

    let mut reversed = coords.clone();
    reversed.reverse();

    let forward = best_rotation(coords);
    let backward = best_rotation(reversed);
    if forward <= backward {
        forward
    } else {
        backward
    }
}

fn canonicalize_polygon(polygon: &Polygon<f64>) -> (Vec<(i64, i64)>, Vec<Vec<(i64, i64)>>) {
    let exterior = canonicalize_ring(polygon.exterior());
    let mut interiors: Vec<_> = polygon.interiors().iter().map(canonicalize_ring).collect();
    interiors.sort();
    (exterior, interiors)
}

fn canonicalize_multipolygon(
    multi: &MultiPolygon<f64>,
) -> Vec<(Vec<(i64, i64)>, Vec<Vec<(i64, i64)>>)> {
    let mut polygons: Vec<_> = multi.0.iter().map(canonicalize_polygon).collect();
    polygons.sort();
    polygons
}

fn assert_polygonize_fixture(name: &str) {
    let input_lines = load_input_lines(name);
    let actual = polygonize(input_lines);
    let expected = load_expected_polygons_from("output", name);

    let actual_canonical = canonicalize_multipolygon(&actual);
    let expected_canonical = canonicalize_multipolygon(&expected);

    assert_eq!(
        actual_canonical, expected_canonical,
        "fixture mismatch: {name}"
    );
}

fn assert_polygonize_fixture_with_output_dir(name: &str, output_kind: &str) {
    let input_lines = load_input_lines(name);
    let actual = polygonize(input_lines);
    let expected = load_expected_polygons_from(output_kind, name);

    let actual_canonical = canonicalize_multipolygon(&actual);
    let expected_canonical = canonicalize_multipolygon(&expected);

    assert_eq!(
        actual_canonical, expected_canonical,
        "fixture mismatch: {name} against {output_kind}"
    );
}

fn assert_polygonize_fixture_with_nodify(name: &str, snap_radius: f64, output_kind: &str) {
    let input_lines = load_input_lines(name);
    let noded_lines = nodify_lines(input_lines.into_iter(), snap_radius);
    let actual = polygonize(noded_lines);
    let expected = load_expected_polygons_from(output_kind, name);

    let actual_canonical = canonicalize_multipolygon(&actual);
    let expected_canonical = canonicalize_multipolygon(&expected);

    assert_eq!(
        actual_canonical, expected_canonical,
        "fixture mismatch: {name} with nodify against {output_kind}"
    );
}

fn assert_polygonize_fixture_with_nodify_exact_region_polygons(
    name: &str,
    snap_radius: f64,
    output_kind: &str,
) {
    let input_lines = load_input_lines(name);
    let noded_lines = nodify_lines(input_lines.into_iter(), snap_radius);
    let actual = polygonize(noded_lines);
    let expected = load_expected_polygons_from(output_kind, name);

    let actual_region_polygons: Vec<_> = canonicalize_multipolygon(&actual)
        .into_iter()
        .filter(|(_exterior, interiors)| interiors.is_empty())
        .collect();
    let expected_region_polygons = canonicalize_multipolygon(&expected);

    assert_eq!(
        actual_region_polygons, expected_region_polygons,
        "fixture mismatch: exact region polygons differ for {name} against {output_kind}"
    );
}

#[test]
fn polygonize_simple_square_builds_single_shell() {
    // Explores the base case: one closed ring should produce exactly one polygon shell.
    assert_polygonize_fixture("simple_square");
}

#[test]
fn polygonize_discards_dangling_spur() {
    // Explores dangle deletion: a spur attached to a valid ring must not create extra polygons.
    assert_polygonize_fixture("square_with_dangle");
}

#[test]
fn polygonize_discards_cut_bridge_between_faces() {
    // Explores cut-edge deletion: a bridge connecting two independent rings should be removed.
    assert_polygonize_fixture("two_squares_with_bridge");
}

#[test]
fn polygonize_assigns_inner_ring_as_hole() {
    // Explores shell-hole assignment: a nested inner ring should become a hole of the outer shell.
    assert_polygonize_fixture("donut_hole");
}

#[test]
fn polygonize_ignores_duplicate_boundary_segments() {
    // Explores duplicate-line robustness: repeated coincident edges should not create extra polygons.
    assert_polygonize_fixture("duplicate_boundary_segments");
}

#[test]
fn polygonize_with_nodify_preserves_pre_noded_fixture_results() {
    // Explores nodify stability: pre-noded inputs should polygonize identically with or without a nodify pre-pass.
    for name in [
        "simple_square",
        "square_with_dangle",
        "two_squares_with_bridge",
        "donut_hole",
        "duplicate_boundary_segments",
    ] {
        assert_polygonize_fixture_with_nodify(name, 1e-9, "output");
    }
}

#[test]
fn polygonize_crossing_lines_requires_nodify_to_form_face() {
    // Explores noding requirement: interior line intersections should only produce a face after nodify splits crossing segments.
    assert_polygonize_fixture_with_output_dir("noding_required_crossing_lines", "output_raw");
    assert_polygonize_fixture_with_nodify("noding_required_crossing_lines", 1e-9, "output_nodify");
}

#[test]
fn polygonize_collinear_overlap_requires_nodify_to_form_face() {
    // Explores collinear overlap noding: partially overlapping collinear edges should only close a face after nodify splits overlap endpoints.
    assert_polygonize_fixture_with_output_dir("collinear_overlap_requires_nodify", "output_raw");
    assert_polygonize_fixture_with_nodify(
        "collinear_overlap_requires_nodify",
        1e-9,
        "output_nodify",
    );
}

#[test]
fn polygonize_multi_level_nesting_retains_nested_hole_hierarchy() {
    // Explores multi-level shell-hole structure: nested rings should produce outer-with-hole, middle-with-hole, and innermost standalone polygon.
    assert_polygonize_fixture("multi_level_nesting");
}

#[test]
fn polygonize_venn_overlaps_split_into_distinct_regions() {
    // Explores overlap partitioning: overlapping polygon inputs should be noded into one polygon per distinct overlap region.
    assert_polygonize_fixture_with_nodify_exact_region_polygons(
        "venn_three_overlapping_rectangles",
        1e-9,
        "output_nodify",
    );
}
