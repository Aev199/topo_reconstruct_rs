//! Development bridge of the editor: serves the built UI (`app/ui/dist`)
//! and forwards `POST /api/<command>` with a JSON body to the same
//! `Service::dispatch` as the Tauri application. Lets the interface be run
//! and tested in a browser without the desktop shell.
//!
//! cargo run --release --example editor_server -- [PORT] [UI_DIR]
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;
use std::path::PathBuf;
use topo_reconstruct_rs::service::Service;

fn main() {
    let mut args = std::env::args().skip(1);
    let port: u16 = args.next().and_then(|p| p.parse().ok()).unwrap_or(8787);
    let ui = PathBuf::from(args.next().unwrap_or_else(|| "app/ui/dist".into()));
    let mut service = Service::new(std::env::temp_dir().join("topo_editor_cache"));
    let listener = TcpListener::bind(("127.0.0.1", port)).expect("bind");
    eprintln!(
        "editor: http://127.0.0.1:{port}/ (UI from {})",
        ui.display()
    );
    for stream in listener.incoming() {
        let Ok(mut stream) = stream else { continue };
        let mut reader = BufReader::new(stream.try_clone().expect("clone"));
        let mut line = String::new();
        if reader.read_line(&mut line).is_err() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let (method, target) = (
            parts.next().unwrap_or(""),
            parts.next().unwrap_or("/").to_string(),
        );
        let mut length = 0usize;
        loop {
            let mut header = String::new();
            if reader.read_line(&mut header).is_err() || header.trim().is_empty() {
                break;
            }
            if let Some((k, v)) = header.split_once(':') {
                if k.eq_ignore_ascii_case("content-length") {
                    length = v.trim().parse().unwrap_or(0);
                }
            }
        }
        let mut body = vec![0; length];
        if reader.read_exact(&mut body).is_err() {
            continue;
        }
        let (status, kind, payload) = if method == "POST" && target.starts_with("/api/") {
            let command = target.trim_start_matches("/api/").to_string();
            let args: serde_json::Value =
                serde_json::from_slice(&body).unwrap_or(serde_json::json!({}));
            let started = std::time::Instant::now();
            let result =
                service.dispatch(&command, args, &mut |stage| eprintln!("  stage {stage}"));
            eprintln!("{command}: {:.2}s", started.elapsed().as_secs_f64());
            let value = match result {
                Ok(v) => serde_json::json!({ "ok": v }),
                Err(e) => serde_json::json!({ "error": e }),
            };
            (
                "200 OK",
                "application/json",
                serde_json::to_vec(&value).unwrap(),
            )
        } else {
            let path = if target == "/" {
                "/index.html".to_string()
            } else {
                target.split('?').next().unwrap().to_string()
            };
            let file = ui.join(path.trim_start_matches('/'));
            let safe = file
                .canonicalize()
                .ok()
                .zip(ui.canonicalize().ok())
                .is_some_and(|(f, root)| f.starts_with(root));
            match std::fs::read(&file) {
                Ok(bytes) if safe => {
                    let kind = match file.extension().and_then(|e| e.to_str()) {
                        Some("html") => "text/html; charset=utf-8",
                        Some("js") => "text/javascript",
                        Some("css") => "text/css",
                        Some("svg") => "image/svg+xml",
                        _ => "application/octet-stream",
                    };
                    ("200 OK", kind, bytes)
                }
                _ => ("404 Not Found", "text/plain", b"not found".to_vec()),
            }
        };
        let _ = write!(
            stream,
            "HTTP/1.1 {status}\r\nContent-Type: {kind}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            payload.len()
        );
        let _ = stream.write_all(&payload);
    }
}
