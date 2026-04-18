use std::fs;
use std::path::PathBuf;

use geo::{Contains, InteriorPoint, Line, LinesIter, MultiPolygon, Polygon, Validation, point};
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

fn ring_has_self_contact(ring: &geo::LineString<f64>, epsilon: f64) -> bool {
    fn orient(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    }

    fn between(value: f64, left: f64, right: f64, epsilon: f64) -> bool {
        let min_value = left.min(right) - epsilon;
        let max_value = left.max(right) + epsilon;
        value >= min_value && value <= max_value
    }

    fn on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64), epsilon: f64) -> bool {
        between(p.0, a.0, b.0, epsilon) && between(p.1, a.1, b.1, epsilon)
    }

    fn segments_contact(
        a1: (f64, f64),
        a2: (f64, f64),
        b1: (f64, f64),
        b2: (f64, f64),
        epsilon: f64,
    ) -> bool {
        let o1 = orient(a1, a2, b1);
        let o2 = orient(a1, a2, b2);
        let o3 = orient(b1, b2, a1);
        let o4 = orient(b1, b2, a2);

        let proper = (o1 > epsilon && o2 < -epsilon || o1 < -epsilon && o2 > epsilon)
            && (o3 > epsilon && o4 < -epsilon || o3 < -epsilon && o4 > epsilon);
        if proper {
            return true;
        }

        (o1.abs() <= epsilon && on_segment(a1, a2, b1, epsilon))
            || (o2.abs() <= epsilon && on_segment(a1, a2, b2, epsilon))
            || (o3.abs() <= epsilon && on_segment(b1, b2, a1, epsilon))
            || (o4.abs() <= epsilon && on_segment(b1, b2, a2, epsilon))
    }

    let mut coordinates: Vec<(f64, f64)> = ring.points().map(|point| (point.x(), point.y())).collect();
    if coordinates.first() == coordinates.last() {
        coordinates.pop();
    }

    let segment_count = coordinates.len();
    if segment_count < 3 {
        return false;
    }

    for first_segment_index in 0..segment_count {
        let first_start = coordinates[first_segment_index];
        let first_end = coordinates[(first_segment_index + 1) % segment_count];

        for second_segment_index in (first_segment_index + 1)..segment_count {
            let first_next_index = (first_segment_index + 1) % segment_count;
            let second_next_index = (second_segment_index + 1) % segment_count;

            let are_adjacent = second_segment_index == first_next_index
                || first_segment_index == second_next_index;
            if are_adjacent {
                continue;
            }

            if first_segment_index == 0 && second_next_index == 0 {
                continue;
            }

            let second_start = coordinates[second_segment_index];
            let second_end = coordinates[second_next_index];

            if segments_contact(first_start, first_end, second_start, second_end, epsilon) {
                return true;
            }
        }
    }

    false
}

fn ring_pair_has_contact(
    ring: &geo::LineString<f64>,
    other_ring: &geo::LineString<f64>,
    epsilon: f64,
) -> bool {
    fn orient(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    }

    fn between(value: f64, left: f64, right: f64, epsilon: f64) -> bool {
        let min_value = left.min(right) - epsilon;
        let max_value = left.max(right) + epsilon;
        value >= min_value && value <= max_value
    }

    fn on_segment(a: (f64, f64), b: (f64, f64), p: (f64, f64), epsilon: f64) -> bool {
        between(p.0, a.0, b.0, epsilon) && between(p.1, a.1, b.1, epsilon)
    }

    fn segments_contact(
        a1: (f64, f64),
        a2: (f64, f64),
        b1: (f64, f64),
        b2: (f64, f64),
        epsilon: f64,
    ) -> bool {
        let o1 = orient(a1, a2, b1);
        let o2 = orient(a1, a2, b2);
        let o3 = orient(b1, b2, a1);
        let o4 = orient(b1, b2, a2);

        let proper = (o1 > epsilon && o2 < -epsilon || o1 < -epsilon && o2 > epsilon)
            && (o3 > epsilon && o4 < -epsilon || o3 < -epsilon && o4 > epsilon);
        if proper {
            return true;
        }

        (o1.abs() <= epsilon && on_segment(a1, a2, b1, epsilon))
            || (o2.abs() <= epsilon && on_segment(a1, a2, b2, epsilon))
            || (o3.abs() <= epsilon && on_segment(b1, b2, a1, epsilon))
            || (o4.abs() <= epsilon && on_segment(b1, b2, a2, epsilon))
    }

    let mut coordinates: Vec<(f64, f64)> = ring.points().map(|point| (point.x(), point.y())).collect();
    let mut other_coordinates: Vec<(f64, f64)> =
        other_ring.points().map(|point| (point.x(), point.y())).collect();

    if coordinates.first() == coordinates.last() {
        coordinates.pop();
    }
    if other_coordinates.first() == other_coordinates.last() {
        other_coordinates.pop();
    }

    if coordinates.len() < 3 || other_coordinates.len() < 3 {
        return false;
    }

    for first_segment_index in 0..coordinates.len() {
        let first_start = coordinates[first_segment_index];
        let first_end = coordinates[(first_segment_index + 1) % coordinates.len()];

        for second_segment_index in 0..other_coordinates.len() {
            let second_start = other_coordinates[second_segment_index];
            let second_end = other_coordinates[(second_segment_index + 1) % other_coordinates.len()];

            if segments_contact(first_start, first_end, second_start, second_end, epsilon) {
                return true;
            }
        }
    }

    false
}

fn polygon_has_pinch_contact(polygon: &Polygon<f64>, epsilon: f64) -> bool {
    if ring_has_self_contact(polygon.exterior(), epsilon) {
        return true;
    }

    if polygon
        .interiors()
        .iter()
        .any(|hole| ring_has_self_contact(hole, epsilon))
    {
        return true;
    }

    if polygon
        .interiors()
        .iter()
        .any(|hole| ring_pair_has_contact(polygon.exterior(), hole, epsilon))
    {
        return true;
    }

    for first_hole_index in 0..polygon.interiors().len() {
        for second_hole_index in (first_hole_index + 1)..polygon.interiors().len() {
            if ring_pair_has_contact(
                &polygon.interiors()[first_hole_index],
                &polygon.interiors()[second_hole_index],
                epsilon,
            ) {
                return true;
            }
        }
    }

    false
}

fn lines_have_positive_collinear_overlap(first: &Line<f64>, second: &Line<f64>, epsilon: f64) -> bool {
    fn cross(a: (f64, f64), b: (f64, f64)) -> f64 {
        a.0 * b.1 - a.1 * b.0
    }

    let first_direction = (first.end.x - first.start.x, first.end.y - first.start.y);
    let second_direction = (second.end.x - second.start.x, second.end.y - second.start.y);

    if first_direction.0.abs() <= epsilon && first_direction.1.abs() <= epsilon {
        return false;
    }

    if cross(first_direction, second_direction).abs() > epsilon {
        return false;
    }

    let second_start_offset = (second.start.x - first.start.x, second.start.y - first.start.y);
    let second_end_offset = (second.end.x - first.start.x, second.end.y - first.start.y);

    if cross(first_direction, second_start_offset).abs() > epsilon
        || cross(first_direction, second_end_offset).abs() > epsilon
    {
        return false;
    }

    let use_x_axis = first_direction.0.abs() >= first_direction.1.abs();
    let (first_min, first_max, second_min, second_max) = if use_x_axis {
        (
            first.start.x.min(first.end.x),
            first.start.x.max(first.end.x),
            second.start.x.min(second.end.x),
            second.start.x.max(second.end.x),
        )
    } else {
        (
            first.start.y.min(first.end.y),
            first.start.y.max(first.end.y),
            second.start.y.min(second.end.y),
            second.start.y.max(second.end.y),
        )
    };

    let overlap_min = first_min.max(second_min);
    let overlap_max = first_max.min(second_max);
    overlap_max - overlap_min > epsilon
}

fn point_is_on_line_segment(point: geo::Coord<f64>, line: &Line<f64>, epsilon: f64) -> bool {
    fn cross(a: (f64, f64), b: (f64, f64)) -> f64 {
        a.0 * b.1 - a.1 * b.0
    }

    let direction = (line.end.x - line.start.x, line.end.y - line.start.y);
    let offset = (point.x - line.start.x, point.y - line.start.y);
    let length_squared = direction.0 * direction.0 + direction.1 * direction.1;

    if length_squared <= epsilon * epsilon {
        return false;
    }

    if cross(direction, offset).abs() > epsilon {
        return false;
    }

    let projection = (offset.0 * direction.0 + offset.1 * direction.1) / length_squared;
    projection >= -epsilon && projection <= 1.0 + epsilon
}

fn line_is_represented_on_output_boundary(
    input_line: &Line<f64>,
    boundary_lines: &[Line<f64>],
    epsilon: f64,
) -> bool {
    if boundary_lines
        .iter()
        .any(|boundary_line| lines_have_positive_collinear_overlap(input_line, boundary_line, epsilon))
    {
        return true;
    }

    let start_on_boundary = boundary_lines
        .iter()
        .any(|boundary_line| point_is_on_line_segment(input_line.start, boundary_line, epsilon));
    let end_on_boundary = boundary_lines
        .iter()
        .any(|boundary_line| point_is_on_line_segment(input_line.end, boundary_line, epsilon));

    start_on_boundary && end_on_boundary
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

#[test]
fn polygonize_failed_hole_probe_point_has_single_owner() {
    let lines = load_input_lines("failed_hole");
    let polygons = polygonize(lines.into_iter());

    let test_point = point! { x:131.85, y: 37.25 };
    let num_polygons_containing_point = polygons.iter().filter(|poly| poly.contains(&test_point)).count();

    assert_eq!(num_polygons_containing_point, 1);
}

#[test]
fn polygonize_touching_hole_overlap_ownership_probe_matches_fixture() {
    assert_polygonize_fixture("touching_hole_overlap_ownership_probe");
}

#[test]
fn polygonize_touching_hole_overlap_ownership_minimal_matches_fixture() {
    assert_polygonize_fixture("touching_hole_overlap_ownership_minimal");
}

#[test]
fn polygonize_touching_hole_overlap_ownership_minimal_invariants() {
    let lines = load_input_lines("touching_hole_overlap_ownership_minimal");
    let polygons = polygonize(lines.clone().into_iter());

    let mut output_boundary_lines: Vec<Line<f64>> = Vec::new();
    for polygon in polygons.iter() {
        output_boundary_lines.extend(polygon.exterior().lines());
        for hole in polygon.interiors() {
            output_boundary_lines.extend(hole.lines());
        }
    }

    for (polygon_index, polygon) in polygons.iter().enumerate() {
        let interior_point = polygon
            .interior_point()
            .unwrap_or_else(|| panic!("polygon {polygon_index} has no interior point"));

        let containing_indices: Vec<_> = polygons
            .iter()
            .enumerate()
            .filter_map(|(candidate_index, candidate)| {
                if candidate.contains(&interior_point) {
                    Some(candidate_index)
                } else {
                    None
                }
            })
            .collect();

        let containing_count = containing_indices.len();

        assert_eq!(
            containing_count, 1,
            "polygon {polygon_index} interior point is contained by {containing_count} polygons"
        );

        assert!(
            !polygon_has_pinch_contact(polygon, 0.0),
            "polygon {polygon_index} has pinch/touch boundary contact"
        );
    }

    for (line_index, input_line) in lines.iter().enumerate() {
        let dx = input_line.end.x - input_line.start.x;
        let dy = input_line.end.y - input_line.start.y;
        if dx.abs() <= 1e-12 && dy.abs() <= 1e-12 {
            continue;
        }

        let appears_on_output_boundary = line_is_represented_on_output_boundary(
            input_line,
            &output_boundary_lines,
            1e-10,
        );

        assert!(
            appears_on_output_boundary,
            "input line {line_index} ({:?} -> {:?}) is not represented on output polygon boundaries",
            input_line.start,
            input_line.end
        );
    }
}

#[test]
fn polygonize_filters_invalid_polygon_from_touching_hole_assignment() {
    fn ring_lines(points: &[(f64, f64)]) -> Vec<Line<f64>> {
        points
            .windows(2)
            .map(|segment| {
                let start = geo::Coord {
                    x: segment[0].0,
                    y: segment[0].1,
                };
                let end = geo::Coord {
                    x: segment[1].0,
                    y: segment[1].1,
                };
                Line::new(start, end)
            })
            .collect()
    }

    let mut lines = Vec::new();
    lines.extend(ring_lines(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ]));
    lines.extend(ring_lines(&[
        (0.0, 4.0),
        (3.0, 4.0),
        (3.0, 6.0),
        (0.0, 6.0),
        (0.0, 4.0),
    ]));

    let polygons = polygonize(lines.into_iter());
    assert!(polygons.0.iter().all(|polygon| polygon.is_valid()));
}

#[test]
fn polygonize_handles_nested_shell_overlap_without_double_containment() {
    let lines = load_input_lines("nested_shell_overlap_minimal");
    let polygons = polygonize(lines.into_iter());
    let test_point = point! { x: 131.85, y: 37.25 };
    let containing_count = polygons
        .0
        .iter()
        .filter(|polygon| polygon.contains(&test_point))
        .count();

    assert_eq!(containing_count, 1);
}
