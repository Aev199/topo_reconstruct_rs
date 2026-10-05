fn main() {
    // A build with GMSH_DLL_PATH (the Gmsh SDK's gmsh-4.x.dll) carries the
    // library inside the executable.
    println!("cargo:rerun-if-env-changed=GMSH_DLL_PATH");
    println!("cargo:rustc-check-cfg=cfg(gmsh_embedded)");
    if std::env::var_os("GMSH_DLL_PATH").is_some_and(|p| std::path::Path::new(&p).is_file()) {
        println!("cargo:rustc-cfg=gmsh_embedded");
    }
    tauri_build::build()
}
