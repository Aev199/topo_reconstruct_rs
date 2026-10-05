// Backend access: the Tauri command `call` in the desktop application, the
// development bridge (`examples/editor_server.rs`) over HTTP in a browser.
// Both forward to the same Rust `Service::dispatch`.
import { invoke } from '@tauri-apps/api/core';
import { listen } from '@tauri-apps/api/event';
import { ask, open, save } from '@tauri-apps/plugin-dialog';

export const desktop = typeof window !== 'undefined' && '__TAURI_INTERNALS__' in window;

export async function call(command, args = {}) {
  if (desktop) {
    return invoke('call', { command, args });
  }
  const response = await fetch(`./api/${command}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(args),
  });
  const body = await response.json();
  if ('error' in body) throw new Error(body.error);
  return body.ok;
}

/// Pipeline stage names while a model opens (desktop only).
export async function onProgress(handler) {
  if (desktop) await listen('progress', (event) => handler(event.payload));
}

const FILTERS = {
  model: [{ name: 'Модель ЛИРА', extensions: ['txt'] }],
  project: [{ name: 'Проект', extensions: ['topo.json', 'json'] }],
  report: [{ name: 'Результат JSON', extensions: ['json'] }],
  plaxis: [{ name: 'Файл обмена PLAXIS', extensions: ['plaxis.json', 'json'] }],
  mxt: [{ name: 'MIDAS Civil', extensions: ['mxt'] }],
};

/// A file path chosen by the user (a dialog on the desktop, a prompt in a
/// browser where the bridge reads local paths).
export async function pickFile(kind, saving = false) {
  if (desktop) {
    const options = { filters: FILTERS[kind] };
    const path = saving ? await save(options) : await open({ ...options, multiple: false });
    return path || null;
  }
  const label = { model: 'Путь к модели ЛИРА', project: 'Путь к проекту', report: 'Путь для результата', plaxis: 'Путь для файла обмена PLAXIS', mxt: 'Путь для файла MIDAS (.mxt)' }[kind];
  return window.prompt(label) || null;
}

/// A yes/no question (a native dialog on the desktop).
export async function question(text) {
  if (desktop) return ask(text, { title: 'Редактор геометрии', kind: 'warning' });
  return window.confirm(text);
}
