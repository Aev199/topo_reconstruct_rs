// Desktop shell of the editor: one command, `call`, forwards to the same
// `Service::dispatch` as the development bridge. Long operations (opening
// a model reconstructs it) run off the UI thread and report their stages
// as `progress` events.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::sync::Mutex;
use tauri::{Emitter, Manager};
use topo_reconstruct_rs::service::Service;

struct Editor(Mutex<Service>);

#[tauri::command]
async fn call(app: tauri::AppHandle, command: String, args: serde_json::Value) -> Result<serde_json::Value, String> {
    tauri::async_runtime::spawn_blocking(move || {
        let editor = app.state::<Editor>();
        let mut service = editor.0.lock().map_err(|_| "the editor state is unavailable".to_string())?;
        service.dispatch(&command, args, &mut |stage| {
            let _ = app.emit("progress", stage);
        })
    })
    .await
    .map_err(|e| e.to_string())?
}

/// Portable: the frame cache and the WebView2 data live in a folder next
/// to the executable when it is writable (else in the user's local app
/// data), so the program runs from any folder or USB drive without setup.
fn data_dir(app: &tauri::App) -> std::path::PathBuf {
    let beside = std::env::current_exe()
        .ok()
        .and_then(|exe| exe.parent().map(|d| d.join("topo-editor-data")));
    if let Some(dir) = beside {
        if std::fs::create_dir_all(&dir).is_ok()
            && std::fs::write(dir.join(".write-test"), b"").is_ok()
        {
            let _ = std::fs::remove_file(dir.join(".write-test"));
            return dir;
        }
    }
    app.path()
        .app_local_data_dir()
        .unwrap_or_else(|_| std::env::temp_dir().join("topo-editor"))
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            let data = data_dir(app);
            app.manage(Editor(Mutex::new(Service::new(data.join("cache")))));
            tauri::WebviewWindowBuilder::new(app, "main", tauri::WebviewUrl::default())
                .title("Topo Editor")
                .inner_size(1400., 900.)
                .min_inner_size(900., 600.)
                .data_directory(data.join("webview"))
                .build()?;
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![call])
        .run(tauri::generate_context!())
        .expect("error while running the editor");
}
