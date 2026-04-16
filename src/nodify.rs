use geo::{
    Coord, GeoFloat, Line,
    algorithm::sweep::Intersections,
    line_intersection::LineIntersection,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug)]
struct TaggedLine<T: GeoFloat> {
    idx: usize,
    line: Line<T>,
}

// `Intersections` needs an input type implementing `Cross`.
// In recent geo versions, the easiest path is usually to wrap your line
// in a small struct like this and implement the required methods.
// The exact trait details may vary slightly by geo version.
impl<T: GeoFloat> geo::algorithm::sweep::Cross for TaggedLine<T> {
    type Scalar = T;

    fn line(&self) -> Line<Self::Scalar> {
        self.line
    }
}

/// Quantized key so we can deduplicate points/subsegments robustly.
/// For production, choose a snap_radius appropriate for your coordinate system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct QCoord(i64, i64);

fn qcoord_raw<T: GeoFloat>(c: Coord<T>, snap_radius: T) -> QCoord {
    QCoord(
        (c.x / snap_radius)
            .round()
            .to_i64()
            .expect("failed to quantize x coordinate"),
        (c.y / snap_radius)
            .round()
            .to_i64()
            .expect("failed to quantize y coordinate"),
    )
}

fn unqcoord_raw<T: GeoFloat>(q: QCoord, snap_radius: T) -> Coord<T> {
    Coord {
        x: T::from(q.0).expect("failed to convert quantized x coordinate") * snap_radius,
        y: T::from(q.1).expect("failed to convert quantized y coordinate") * snap_radius,
    }
}

struct Quantizer<T: GeoFloat> {
    snap_radius: T,
    first_seen: BTreeMap<QCoord, Coord<T>>,
}

impl<T: GeoFloat> Quantizer<T> {
    fn new(snap_radius: T) -> Self {
        Self {
            snap_radius,
            first_seen: BTreeMap::new(),
        }
    }

    fn qcoord(&mut self, c: Coord<T>) -> QCoord {
        let q = qcoord_raw(c, self.snap_radius);
        self.first_seen.entry(q).or_insert(c);
        q
    }

    fn unqcoord(&self, q: QCoord) -> Coord<T> {
        self.first_seen
            .get(&q)
            .copied()
            .unwrap_or_else(|| unqcoord_raw(q, self.snap_radius))
    }
}

fn normalize_segment(qa: QCoord, qb: QCoord) -> (QCoord, QCoord) {
    if qa <= qb { (qa, qb) } else { (qb, qa) }
}

pub fn nodify_lines<T: GeoFloat>(
    input: impl IntoIterator<Item = Line<T>>,
    snap_radius: T,
) -> Vec<Line<T>> {
    assert!(snap_radius > T::zero(), "snap_radius must be > 0");

    let mut quantizer = Quantizer::new(snap_radius);

    let tagged: Vec<TaggedLine<T>> = input
        .into_iter()
        .enumerate()
        .map(|(idx, line)| TaggedLine { idx, line })
        .collect();

    // Per original line idx, the set of cut points that lie on that line.
    let mut cuts: BTreeMap<usize, BTreeSet<QCoord>> = BTreeMap::new();

    // Seed each line with its endpoints.
    for tl in &tagged {
        cuts.entry(tl.idx)
            .or_default()
            .insert(quantizer.qcoord(tl.line.start));
        cuts.entry(tl.idx)
            .or_default()
            .insert(quantizer.qcoord(tl.line.end));
    }

    // Discover all pairwise interactions with the sweep.
    for (a, b, ix) in Intersections::from_iter(tagged.iter().copied()) {
        match ix {
            LineIntersection::SinglePoint { intersection, .. } => {
                let q = quantizer.qcoord(intersection);
                cuts.entry(a.idx).or_default().insert(q);
                cuts.entry(b.idx).or_default().insert(q);
            }
            LineIntersection::Collinear { intersection } => {
                // Overlap segment endpoints become cut points on both originals.
                let q_start = quantizer.qcoord(intersection.start);
                let q_end = quantizer.qcoord(intersection.end);
                cuts.entry(a.idx).or_default().insert(q_start);
                cuts.entry(a.idx).or_default().insert(q_end);
                cuts.entry(b.idx).or_default().insert(q_start);
                cuts.entry(b.idx).or_default().insert(q_end);
            }
        }
    }

    // Split each original line at its ordered cut points.
    let mut out: BTreeSet<(QCoord, QCoord)> = BTreeSet::new();

    for tl in &tagged {
        let Some(pts) = cuts.get(&tl.idx) else { continue };

        let mut pts: Vec<Coord<T>> = pts.iter().map(|&q| quantizer.unqcoord(q)).collect();

        let dx = tl.line.end.x - tl.line.start.x;
        let dy = tl.line.end.y - tl.line.start.y;
        pts.sort_by(|p1, p2| {
            let t1 = (p1.x - tl.line.start.x) * dx + (p1.y - tl.line.start.y) * dy;
            let t2 = (p2.x - tl.line.start.x) * dx + (p2.y - tl.line.start.y) * dy;
            t1.total_cmp(&t2)
        });

        pts.dedup_by(|a, b| qcoord_raw(*a, snap_radius) == qcoord_raw(*b, snap_radius));

        for w in pts.windows(2) {
            let p0 = w[0];
            let p1 = w[1];
            let q0 = qcoord_raw(p0, snap_radius);
            let q1 = qcoord_raw(p1, snap_radius);

            // Skip zero-length pieces
            if q0 == q1 {
                continue;
            }

            out.insert(normalize_segment(q0, q1));
        }
    }

    let lines: Vec<_> = out
        .into_iter()
        .map(|(a, b)| Line::new(quantizer.unqcoord(a), quantizer.unqcoord(b)))
        .collect();

    lines
}