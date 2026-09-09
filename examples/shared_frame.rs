//! Run: cargo run --example shared_frame > frame.json
//! Kernel demonstration only: no FE recognition or meshing is claimed.
use topo_reconstruct_rs::reconstruction::{Model, PlaneFrame};
fn main() {
    let mut model = Model::new(1e-8, 0.03).unwrap();
    let slab = model.add_plane(PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap());
    let wall = model.add_plane(PlaneFrame::new([0.; 3], [0., 1., 0.]).unwrap());
    let ids = [
        [0., 0., 0.],
        [4., 0., 0.],
        [4., 3., 0.],
        [0., 3., 0.],
        [0., 0., 3.],
        [4., 0., 3.],
    ]
    .map(|p| model.add_vertex(p).unwrap());
    model
        .add_surface(slab, vec![vec![ids[0], ids[1], ids[2], ids[3]]], vec![1])
        .unwrap();
    model
        .add_surface(wall, vec![vec![ids[1], ids[0], ids[4], ids[5]]], vec![2])
        .unwrap();
    println!("{}", serde_json::to_string_pretty(&model).unwrap());
}
