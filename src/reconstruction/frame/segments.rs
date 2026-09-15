//! Constructive spans, independent of source FE subdivision and property changes.
use super::*;

pub(super) fn split(
    mesh: &MeshData,
    report: &recognize::Report,
    planes: &planes::Report,
) -> Vec<(recognize::RecognizedAxis, bool)> {
    let mut incidence = BTreeMap::<u32, BTreeSet<usize>>::new();
    for (i, axis) in report.axes.iter().enumerate() {
        for a in &axis.anchors {
            incidence.entry(a.node).or_default().insert(i);
        }
    }
    let mut owners = BTreeMap::<u32, BTreeSet<usize>>::new();
    for (i, plane) in planes.patches.iter().enumerate() {
        for &n in &plane.source_nodes {
            owners.entry(n).or_default().insert(i);
        }
    }
    let mut result = vec![];
    for axis in &report.axes {
        let mut anchors = axis.anchors.clone();
        anchors.sort_by(|a, b| a.t.total_cmp(&b.t));
        let direction =
            (mesh.nodes[&axis.endpoint_nodes[1]] - mesh.nodes[&axis.endpoint_nodes[0]]).normalize();
        let mut cuts = vec![0];
        for i in 1..anchors.len() - 1 {
            let n = anchors[i].node;
            let contact = incidence[&n].len() > 1
                || owners.get(&n).is_some_and(|ps| {
                    ps.iter().any(|&p| {
                        let transverse = DVec3::from_array(planes.patches[p].plane.normal)
                            .dot(direction)
                            .abs()
                            > report.policy.angle.sin();
                        // A continuous coplanar connection has only entry/exit cuts.
                        transverse
                            || !owners
                                .get(&anchors[i - 1].node)
                                .is_some_and(|s| s.contains(&p))
                            || !owners
                                .get(&anchors[i + 1].node)
                                .is_some_and(|s| s.contains(&p))
                    })
                });
            if contact {
                cuts.push(i);
            }
        }
        cuts.push(anchors.len() - 1);
        let segmented = cuts.len() > 2;
        for w in cuts.windows(2) {
            let lo = anchors[w[0]].t;
            let hi = anchors[w[1]].t;
            let mut part = axis.clone();
            part.endpoint_nodes = [anchors[w[0]].node, anchors[w[1]].node];
            part.endpoints = part.endpoint_nodes.map(|n| mesh.nodes[&n].to_array());
            part.anchors = anchors[w[0]..=w[1]].to_vec();
            for a in &mut part.anchors {
                a.t = (a.t - lo) / (hi - lo);
            }
            part.spans = axis
                .spans
                .iter()
                .filter(|s| s.start_t >= lo && s.end_t <= hi)
                .cloned()
                .map(|mut s| {
                    s.start_t = (s.start_t - lo) / (hi - lo);
                    s.end_t = (s.end_t - lo) / (hi - lo);
                    s
                })
                .collect();
            result.push((part, segmented));
        }
    }
    result
}
