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

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            let cache = app.path().app_cache_dir()?;
            app.manage(Editor(Mutex::new(Service::new(cache))));
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![call])
        .run(tauri::generate_context!())
        .expect("error while running the editor");
}
