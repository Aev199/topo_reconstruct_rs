//! Gmsh through its C API, loaded at run time (`libgmsh.so` / `gmsh-4.x.dll`):
//! the application meshes planar surfaces with shared boundary curves, so
//! neighbouring surfaces get conforming meshes, and lines (bars) embedded in
//! them or standing alone.
//!
//! The library is searched in `TOPO_GMSH_LIB`, beside the executable and in
//! the working directory. Gmsh keeps global state, so meshing is serialized.
use libloading::{Library, Symbol};
use std::ffi::{c_char, c_double, c_int, c_void, CString};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

static LOCK: Mutex<()> = Mutex::new(());

#[derive(Debug, Clone, PartialEq)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl std::error::Error for Error {}

/// A planar surface: loops of signed line numbers (1-based into
/// `Input::lines`, negative = traversed backwards), the first loop outer, the
/// others holes; lines lying inside the surface that its mesh must follow.
#[derive(Debug, Clone, Default)]
pub struct SurfaceInput {
    pub loops: Vec<Vec<i64>>,
    pub embedded_lines: Vec<usize>,
}

#[derive(Debug, Clone, Default)]
pub struct Input {
    pub vertices: Vec<[f64; 3]>,
    /// Straight lines between vertex indices.
    pub lines: Vec<[usize; 2]>,
    pub surfaces: Vec<SurfaceInput>,
    /// Vertices that lie inside a surface (point loads, bar ends).
    pub embedded_points: Vec<(usize, usize)>,
    /// Target element size.
    pub size: f64,
    /// Quadrilaterals where possible (recombination).
    pub quads: bool,
}

#[derive(Debug, Clone, Default)]
pub struct Output {
    pub nodes: Vec<[f64; 3]>,
    /// Node of every input vertex.
    pub vertex_nodes: Vec<usize>,
    /// Per input surface: elements as 3 or 4 node indices.
    pub surface_elements: Vec<Vec<Vec<usize>>>,
    /// Per input line: its segments as node pairs, in the line's direction.
    pub line_elements: Vec<Vec<[usize; 2]>>,
}

type Ierr = *mut c_int;

pub struct Gmsh {
    lib: Library,
    path: PathBuf,
}

/// Where the library may be.
fn candidates() -> Vec<PathBuf> {
    let mut found = vec![];
    if let Some(p) = std::env::var_os("TOPO_GMSH_LIB") {
        found.push(PathBuf::from(p));
    }
    let mut dirs = vec![];
    if let Some(dir) = std::env::current_exe().ok().and_then(|e| e.parent().map(Path::to_path_buf)) {
        dirs.push(dir);
    }
    if let Ok(dir) = std::env::current_dir() {
        dirs.push(dir);
    }
    for dir in dirs {
        if let Ok(entries) = std::fs::read_dir(&dir) {
            let mut names: Vec<PathBuf> = entries
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_lowercase();
                    (name.starts_with("gmsh") && name.ends_with(".dll")) || name.starts_with("libgmsh.so") || name.starts_with("libgmsh.dylib")
                })
                .collect();
            names.sort();
            found.extend(names);
        }
    }
    found
}

macro_rules! call {
    ($lib:expr, $name:literal, fn($($t:ty),*) $(-> $r:ty)?, $($a:expr),*) => {{
        let f: Symbol<unsafe extern "C" fn($($t),*) $(-> $r)?> =
            $lib.get(concat!($name, "\0").as_bytes()).map_err(|e| Error(format!("{}: {e}", $name)))?;
        f($($a),*)
    }};
}

impl Gmsh {
    /// The library, or why it cannot be used.
    pub fn load() -> Result<Gmsh, Error> {
        let mut last = Error("Gmsh library not found: put gmsh-4.x.dll (libgmsh.so) beside the program or set TOPO_GMSH_LIB".into());
        for path in candidates() {
            // SAFETY: loading a shared library runs its initializers; it is Gmsh's own.
            match unsafe { Library::new(&path) } {
                Ok(lib) => return Ok(Gmsh { lib, path }),
                Err(e) => last = Error(format!("{}: {e}", path.display())),
            }
        }
        Err(last)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    fn check(&self, ierr: c_int, what: &str) -> Result<(), Error> {
        if ierr != 0 {
            Err(Error(format!("Gmsh: {what} failed ({ierr})")))
        } else {
            Ok(())
        }
    }

    /// Mesh the input.
    pub fn mesh(&self, input: &Input) -> Result<Output, Error> {
        let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: calls into the Gmsh C API with pointers valid for each call;
        // arrays returned by Gmsh are copied and freed with `gmshFree`.
        unsafe { self.mesh_locked(input) }
    }

    unsafe fn mesh_locked(&self, input: &Input) -> Result<Output, Error> {
        let lib = &self.lib;
        let mut ierr: c_int = 0;
        call!(lib, "gmshInitialize", fn(c_int, *mut *mut c_char, c_int, c_int, Ierr), 0, std::ptr::null_mut(), 0, 0, &mut ierr);
        self.check(ierr, "initialize")?;
        let result = self.mesh_initialized(input);
        let mut e2: c_int = 0;
        call!(lib, "gmshFinalize", fn(Ierr), &mut e2);
        result
    }

    unsafe fn mesh_initialized(&self, input: &Input) -> Result<Output, Error> {
        let lib = &self.lib;
        let mut ierr: c_int = 0;
        let option = |name: &str, value: f64| -> Result<(), Error> {
            let mut e: c_int = 0;
            let c = CString::new(name).unwrap();
            call!(lib, "gmshOptionSetNumber", fn(*const c_char, c_double, Ierr), c.as_ptr(), value, &mut e);
            self.check(e, name)
        };
        option("General.Terminal", 0.)?;
        option("Mesh.MeshSizeMax", input.size)?;
        option("Mesh.MeshSizeMin", input.size * 0.2)?;
        // Frontal-Delaunay for triangles; quads by recombination.
        option("Mesh.Algorithm", 6.)?;
        if input.quads {
            option("Mesh.RecombineAll", 1.)?;
            option("Mesh.RecombinationAlgorithm", 2.)?;
            option("Mesh.Algorithm", 8.)?;
        }
        let name = CString::new("topo").unwrap();
        call!(lib, "gmshModelAdd", fn(*const c_char, Ierr), name.as_ptr(), &mut ierr);
        self.check(ierr, "model")?;
        for (i, v) in input.vertices.iter().enumerate() {
            call!(lib, "gmshModelGeoAddPoint", fn(c_double, c_double, c_double, c_double, c_int, Ierr) -> c_int,
                v[0], v[1], v[2], input.size, (i + 1) as c_int, &mut ierr);
            self.check(ierr, "point")?;
        }
        for (i, l) in input.lines.iter().enumerate() {
            call!(lib, "gmshModelGeoAddLine", fn(c_int, c_int, c_int, Ierr) -> c_int,
                (l[0] + 1) as c_int, (l[1] + 1) as c_int, (i + 1) as c_int, &mut ierr);
            self.check(ierr, "line")?;
        }
        let mut surface_tag = 0;
        let mut loop_tag = 0;
        let mut surface_tags = vec![];
        for s in &input.surfaces {
            let mut wires = vec![];
            for l in &s.loops {
                loop_tag += 1;
                let curves: Vec<c_int> = l.iter().map(|&c| c as c_int).collect();
                call!(lib, "gmshModelGeoAddCurveLoop", fn(*const c_int, usize, c_int, c_int, Ierr) -> c_int,
                    curves.as_ptr(), curves.len(), loop_tag, 0, &mut ierr);
                self.check(ierr, "curve loop")?;
                wires.push(loop_tag);
            }
            surface_tag += 1;
            call!(lib, "gmshModelGeoAddPlaneSurface", fn(*const c_int, usize, c_int, Ierr) -> c_int,
                wires.as_ptr(), wires.len(), surface_tag, &mut ierr);
            self.check(ierr, "plane surface")?;
            surface_tags.push(surface_tag);
        }
        call!(lib, "gmshModelGeoSynchronize", fn(Ierr), &mut ierr);
        self.check(ierr, "synchronize")?;
        for (s, surface) in input.surfaces.iter().enumerate() {
            if !surface.embedded_lines.is_empty() {
                let tags: Vec<c_int> = surface.embedded_lines.iter().map(|&l| (l + 1) as c_int).collect();
                call!(lib, "gmshModelMeshEmbed", fn(c_int, *const c_int, usize, c_int, c_int, Ierr),
                    1, tags.as_ptr(), tags.len(), 2, surface_tags[s], &mut ierr);
                self.check(ierr, "embed line")?;
            }
        }
        for &(vertex, s) in &input.embedded_points {
            let tag = (vertex + 1) as c_int;
            call!(lib, "gmshModelMeshEmbed", fn(c_int, *const c_int, usize, c_int, c_int, Ierr),
                0, &tag, 1, 2, surface_tags[s], &mut ierr);
            self.check(ierr, "embed point")?;
        }
        let dim = if input.surfaces.is_empty() { 1 } else { 2 };
        call!(lib, "gmshModelMeshGenerate", fn(c_int, Ierr), dim, &mut ierr);
        self.check(ierr, "generate")?;

        // All nodes.
        let (tags, coords) = self.nodes(-1, -1)?;
        let mut index: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        let mut nodes = Vec::with_capacity(tags.len());
        for (i, &t) in tags.iter().enumerate() {
            index.insert(t, i);
            nodes.push([coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]);
        }
        let mut out = Output { nodes, ..Default::default() };
        for i in 0..input.vertices.len() {
            let (t, _) = self.nodes(0, (i + 1) as c_int)?;
            let node = t.first().and_then(|t| index.get(t)).copied().ok_or_else(|| Error(format!("Gmsh: vertex {i} has no node")))?;
            out.vertex_nodes.push(node);
        }
        for l in 0..input.lines.len() {
            let mut segments = vec![];
            for (kind, flat) in self.elements(1, (l + 1) as c_int)? {
                if kind == 1 {
                    for pair in flat.chunks(2) {
                        segments.push([index[&pair[0]], index[&pair[1]]]);
                    }
                }
            }
            out.line_elements.push(segments);
        }
        for &tag in &surface_tags {
            let mut elements = vec![];
            for (kind, flat) in self.elements(2, tag)? {
                let n = match kind {
                    2 => 3,
                    3 => 4,
                    _ => continue,
                };
                for e in flat.chunks(n) {
                    elements.push(e.iter().map(|t| index[t]).collect());
                }
            }
            out.surface_elements.push(elements);
        }
        Ok(out)
    }

    /// Node tags and coordinates of an entity (-1, -1: the whole model).
    unsafe fn nodes(&self, dim: c_int, tag: c_int) -> Result<(Vec<usize>, Vec<f64>), Error> {
        let lib = &self.lib;
        let (mut tags, mut tags_n): (*mut usize, usize) = (std::ptr::null_mut(), 0);
        let (mut coord, mut coord_n): (*mut f64, usize) = (std::ptr::null_mut(), 0);
        let (mut par, mut par_n): (*mut f64, usize) = (std::ptr::null_mut(), 0);
        let mut ierr: c_int = 0;
        call!(lib, "gmshModelMeshGetNodes", fn(*mut *mut usize, *mut usize, *mut *mut f64, *mut usize, *mut *mut f64, *mut usize, c_int, c_int, c_int, c_int, Ierr),
            &mut tags, &mut tags_n, &mut coord, &mut coord_n, &mut par, &mut par_n, dim, tag, 1, 0, &mut ierr);
        self.check(ierr, "get nodes")?;
        let out = (
            std::slice::from_raw_parts(tags, tags_n).to_vec(),
            std::slice::from_raw_parts(coord, coord_n).to_vec(),
        );
        self.free(tags as *mut c_void)?;
        self.free(coord as *mut c_void)?;
        self.free(par as *mut c_void)?;
        Ok(out)
    }

    /// Elements of an entity: (gmsh type, node tags of all its elements).
    unsafe fn elements(&self, dim: c_int, tag: c_int) -> Result<Vec<(c_int, Vec<usize>)>, Error> {
        let lib = &self.lib;
        let (mut types, mut types_n): (*mut c_int, usize) = (std::ptr::null_mut(), 0);
        let (mut etags, mut etags_n, mut etags_nn): (*mut *mut usize, *mut usize, usize) = (std::ptr::null_mut(), std::ptr::null_mut(), 0);
        let (mut ntags, mut ntags_n, mut ntags_nn): (*mut *mut usize, *mut usize, usize) = (std::ptr::null_mut(), std::ptr::null_mut(), 0);
        let mut ierr: c_int = 0;
        call!(lib, "gmshModelMeshGetElements",
            fn(*mut *mut c_int, *mut usize, *mut *mut *mut usize, *mut *mut usize, *mut usize, *mut *mut *mut usize, *mut *mut usize, *mut usize, c_int, c_int, Ierr),
            &mut types, &mut types_n, &mut etags, &mut etags_n, &mut etags_nn, &mut ntags, &mut ntags_n, &mut ntags_nn, dim, tag, &mut ierr);
        self.check(ierr, "get elements")?;
        let mut out = vec![];
        for i in 0..types_n {
            let kind = *types.add(i);
            let n = *ntags_n.add(i);
            out.push((kind, std::slice::from_raw_parts(*ntags.add(i), n).to_vec()));
            self.free(*ntags.add(i) as *mut c_void)?;
            self.free(*etags.add(i) as *mut c_void)?;
        }
        self.free(types as *mut c_void)?;
        self.free(etags as *mut c_void)?;
        self.free(etags_n as *mut c_void)?;
        self.free(ntags as *mut c_void)?;
        self.free(ntags_n as *mut c_void)?;
        Ok(out)
    }

    unsafe fn free(&self, p: *mut c_void) -> Result<(), Error> {
        if !p.is_null() {
            call!(&self.lib, "gmshFree", fn(*mut c_void), p);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A square with a square hole, a neighbouring square sharing one side
    /// and a bar along the shared side: the meshes conform.
    #[test]
    fn conforming_mesh_of_two_surfaces_with_a_hole() {
        let gmsh = match Gmsh::load() {
            Ok(g) => g,
            Err(e) => {
                eprintln!("Gmsh library not available ({e}): test skipped");
                return;
            }
        };
        let v = |x: f64, y: f64| [x, y, 0.];
        let input = Input {
            vertices: vec![
                v(0., 0.), v(4., 0.), v(4., 4.), v(0., 4.), // 0..3 outer
                v(1., 1.), v(3., 1.), v(3., 3.), v(1., 3.), // 4..7 hole
                v(8., 0.), v(8., 4.),                      // 8, 9 neighbour
            ],
            lines: vec![
                [0, 1], [1, 2], [2, 3], [3, 0],     // 1..4
                [4, 5], [5, 6], [6, 7], [7, 4],     // 5..8
                [1, 8], [8, 9], [9, 2],             // 9..11
            ],
            surfaces: vec![
                SurfaceInput { loops: vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]], embedded_lines: vec![] },
                SurfaceInput { loops: vec![vec![9, 10, 11, -2]], embedded_lines: vec![] },
            ],
            embedded_points: vec![],
            size: 1.,
            quads: false,
        };
        let out = gmsh.mesh(&input).expect("mesh");
        assert_eq!(out.surface_elements.len(), 2);
        assert!(out.surface_elements.iter().all(|s| s.len() > 8));
        // The shared line 2 has the same segments from both sides: every node
        // of it is a node of both surfaces.
        let on = |s: usize| -> std::collections::BTreeSet<usize> { out.surface_elements[s].iter().flatten().copied().collect() };
        let shared: std::collections::BTreeSet<usize> = out.line_elements[1].iter().flatten().copied().collect();
        assert!(shared.len() >= 4);
        assert!(shared.iter().all(|n| on(0).contains(n) && on(1).contains(n)));
        // No element in the hole.
        for e in &out.surface_elements[0] {
            let c = e.iter().map(|&n| out.nodes[n]).fold([0.; 2], |a, p| [a[0] + p[0], a[1] + p[1]]);
            let (cx, cy) = (c[0] / e.len() as f64, c[1] / e.len() as f64);
            assert!(!(cx > 1.05 && cx < 2.95 && cy > 1.05 && cy < 2.95), "element in the hole at {cx} {cy}");
        }
        // Total area = 16 - 4 + 16.
        let area: f64 = out.surface_elements.iter().flatten().map(|e| {
            let p: Vec<[f64; 3]> = e.iter().map(|&n| out.nodes[n]).collect();
            let mut a = 0.;
            for i in 1..p.len() - 1 {
                a += ((p[i][0] - p[0][0]) * (p[i + 1][1] - p[0][1]) - (p[i][1] - p[0][1]) * (p[i + 1][0] - p[0][0])).abs() / 2.;
            }
            a
        }).sum();
        assert!((area - 28.).abs() < 1e-6, "{area}");
    }
}
