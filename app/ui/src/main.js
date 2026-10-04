// Editor shell: toolbar, audit findings, selection and edit actions. All
// geometry and every check live in the Rust core; this file only shows
// them and sends edits.
import { call, onProgress, pickFile, question } from './api.js';
import { Viewer } from './viewer.js';

const $ = (id) => document.getElementById(id);
const KIND = {
  invalid_surface: 'Неверный контур поверхности',
  coplanar_overlap: 'Наложение поверхностей в одной плоскости',
  unrepresented_crossing: 'Пересечение поверхностей без общего ребра',
  unrepresented_t_junction: 'Т-примыкание без общего ребра',
  unrepresented_boundary_junction: 'Общая граница без общего ребра',
  unshared_point_contact: 'Точка касания без общей вершины',
  overlapping_bars: 'Наложение стержней',
  unshared_bar_intersection: 'Пересечение стержней без общего узла',
  unshared_bar_surface_intersection: 'Стержень пересекает поверхность без узла',
  bar_in_surface_without_contact: 'Стержень в поверхности без связи',
  gap: 'Зазор',
  short_edge: 'Короткое ребро',
  short_bar_piece: 'Короткий участок стержня',
  sharp_corner: 'Острый угол',
  narrow_face: 'Узкое место поверхности',
  accepted_joint: 'Шов (принят)',
  surface_near_miss: 'Близко к поверхности',
  near_parallel_faces: 'Близкие параллельные грани',
  bar_near_miss: 'Стержни рядом',
  bar_surface_near_miss: 'Стержень рядом с поверхностью',
  short_bar: 'Короткий стержень',
};
const CLASS = { failure: 'ошибка', plaxis: 'PLAXIS', review: 'обзор' };
const OP = {
  move_vertex: 'Сдвиг вершины',
  merge_vertices: 'Слияние вершин',
  delete_surface: 'Удаление поверхности',
  join_surfaces: 'Объединение поверхностей',
  split_edge: 'Разбиение ребра',
  close_gap: 'Закрытие зазора',
  mark_joint: 'Зазор принят как шов',
};
/// Russian text of an edit outcome reported by the core.
function outcome(text) {
  if (!text) return 'готово';
  const rules = [
    [/^moved ([0-9.]+)$/, (m) => `сдвиг ${m[1]} м`],
    [/^merged, moved ([0-9.]+)$/, (m) => `слито, сдвиг ${m[1]} м`],
    [/^deleted$/, () => 'удалено'],
    [/^joined$/, () => 'объединено'],
    [/^split at ([0-9.]+), vertex (\d+)$/, (m) => `разбито в точке ${m[1]} длины, вершина ${m[2]}`],
    [/^marked as joint$/, () => 'принят как шов'],
    [/^(\w+), moved ([0-9.]+)$/, (m) => `${m[1] === 'onto_edge' ? 'на ребро' : m[1] === 'onto_vertex' ? 'в вершину' : 'на поверхность'}, сдвиг ${m[2]} м`],
  ];
  for (const [re, f] of rules) { const m = text.match(re); if (m) return f(m); }
  return text;
}

const state = {
  summary: null,
  scene: null,
  audit: null,
  selection: null, // { kind, ... }
  pending: null, // a second pick an action waits for
  finding: null,
};

const viewer = new Viewer($('view'), onPick);

/// Readable reasons of refused edits (the core reports them as codes).
function explain(message) {
  const m = String(message);
  const move = m.match(/movement_beyond_tolerance: ([0-9.e+-]+)/);
  if (move) return `Отклонено: вершину пришлось бы сдвинуть на ${format(Number(move[1]), 3)} м — больше зазора. Сдвиньте или слейте вершину вручную.`;
  const table = [
    ['bend_InvalidRing', 'Отклонено: закрытие ломает контур поверхности.'],
    ['NonPlanar', 'Отклонено: вершина ушла бы с плоскости одной из своих поверхностей.'],
    ['ShortEdge', 'Отклонено: появилось бы ребро короче 1 мм.'],
    ['InvalidRing', 'Отклонено: контур стал бы неверным (самопересечение или разрыв).'],
    ['interior_bar_anchor', 'Отклонено: это промежуточный узел стержня — стержень согнулся бы.'],
    ['shared_bar_anchor', 'Отклонено: узел общий для нескольких стержней.'],
    ['different stiffness', 'Отклонено: у поверхностей разная жёсткость — материалы не объединяются.'],
    ['no_gap', 'Зазора больше нет.'],
    ['inconsistent_planes', 'Отклонено: плоскости вершины не пересекаются в одной точке.'],
  ];
  for (const [code, text] of table) if (m.includes(code)) return `${text} (${m})`;
  return m;
}

function status(text, error = false) {
  $('status').textContent = text;
  $('status').classList.toggle('error', error);
}

async function busy(text, work) {
  $('busy-text').textContent = text;
  $('busy').hidden = false;
  try {
    return await work();
  } catch (e) {
    status(explain(e.message || e), true);
    return null;
  } finally {
    $('busy').hidden = true;
  }
}

onProgress((stage) => { $('busy-text').textContent = `Восстановление: ${stage}…`; });

function format(x, digits = 3) {
  return Number(x).toFixed(digits);
}

function describe(f) {
  const unit = f.kind === 'sharp_corner' ? '°' : f.kind.includes('overlap') && f.kind.startsWith('coplanar') ? ' м²' : ' м';
  const what = [];
  if ((f.kind === 'gap' || f.kind.endsWith('near_miss') || f.kind === 'accepted_joint' || f.kind === 'unshared_point_contact')
      && f.vertex !== undefined && f.surfaces.length > 1) {
    // A vertex of some surfaces near (or on) another surface.
    const target = f.surfaces[f.surfaces.length - 1];
    return `верш. ${f.vertex} (пов. ${f.surfaces.slice(0, -1).join(', ')}) → пов. ${target}`
      + `${f.value ? ` · ${format(f.value, 4)} м` : ''}`;
  }
  if (f.surfaces.length) what.push(`пов. ${f.surfaces.join(', ')}`);
  if (f.bars.length) what.push(`ст. ${f.bars.join(', ')}`);
  if (f.vertex !== undefined) what.push(`верш. ${f.vertex}`);
  return `${what.join(' · ')}${f.value ? ` · ${format(f.value, 4)}${unit}` : ''}`;
}

async function refresh(scene = true) {
  state.summary = await call('summary');
  state.audit = await call('audit');
  if (scene) {
    state.scene = await call('scene');
    viewer.setScene(state.scene);
  }
  viewer.setFindings(visibleFindings());
  renderSummary();
  renderFindings();
  renderJournal();
  select(null);
}

function visibleFindings() {
  if (!state.audit) return [];
  const show = { failure: $('show-failure').checked, plaxis: $('show-plaxis').checked, review: $('show-review').checked };
  return state.audit.findings.filter((f) => show[f.class]);
}

function renderSummary() {
  const s = state.summary;
  const a = s.audit;
  const accepted = a.counts?.accepted_joint || 0;
  const box = $('audit-summary');
  box.className = !a.passed ? 'failed' : a.plaxis_passed ? 'passed' : 'warning';
  // Two separate verdicts: valid connected geometry, then the PLAXIS
  // profile of small features (accepted joints are exceptions, not items).
  box.innerHTML = `<b>Геометрия и связность: ${a.passed ? 'без ошибок' : `ошибок ${a.failures}`}</b><br>`
    + `<b>Профиль PLAXIS: ${!a.passed ? 'не проверяется до исправления ошибок' : a.plaxis_passed ? 'пройден' : `замечаний ${a.plaxis}`}</b><br>`
    + `${accepted ? `принятых швов ${accepted} · ` : ''}обзор ${a.review - accepted}<br>`
    + `<span class="meta">${s.surfaces} поверхностей · ${s.bars} стержней · правок ${s.edits}${s.dirty ? ' (не сохранены)' : ''}</span>`;
  $('save-project').disabled = false;
  $('export-report').disabled = false;
  $('undo').disabled = s.edits === 0;
  $('redo').disabled = !s.can_redo;
}

function renderFindings() {
  const list = $('findings');
  list.innerHTML = '';
  for (const f of visibleFindings()) {
    const li = document.createElement('li');
    li.className = f.class;
    li.innerHTML = `${KIND[f.kind] || f.kind} <span class="meta">(${CLASS[f.class]})</span><div class="meta">${describe(f)}</div>`;
    li.onclick = () => selectFinding(f, li);
    list.appendChild(li);
  }
}

function renderJournal() {
  const list = $('journal');
  list.innerHTML = '';
  // The journal is fetched lazily: the summary carries its length.
  call('journal').then((entries) => {
    for (const e of entries) {
      const li = document.createElement('li');
      li.textContent = `${OP[e.edit.op] || e.edit.op} — ${outcome(e.outcome)}`;
      if (e.note) { const n = document.createElement('div'); n.className = 'note'; n.textContent = KIND[e.note] || e.note; li.appendChild(n); }
      list.appendChild(li);
    }
  });
}

function selectFinding(f, li) {
  for (const x of document.querySelectorAll('#findings li.selected')) x.classList.remove('selected');
  li?.classList.add('selected');
  state.finding = f;
  const extent = f.points.length > 1 ? Math.hypot(...f.points[0].map((x, k) => x - f.points[1][k])) : f.value;
  viewer.flyTo(f.points[0], extent);
  select({ kind: 'finding', finding: f });
}

function select(sel) {
  state.selection = sel;
  const view = { surfaces: [], vertices: [], edges: [], bars: [] };
  let text = 'ничего';
  if (sel?.kind === 'surface') {
    const s = state.scene.surfaces[sel.surface];
    view.surfaces.push(sel.surface);
    text = `Поверхность ${sel.surface}\nжёсткость ${s.stiffness}, площадь ${format(s.area, 3)} м², КЭ ${s.source_elements}`;
  } else if (sel?.kind === 'vertex') {
    view.vertices.push(sel.vertex);
    const p = state.scene.vertices[sel.vertex];
    text = `Вершина ${sel.vertex}\n${p.map((x) => format(x, 4)).join(', ')}`;
  } else if (sel?.kind === 'edge') {
    view.edges.push(sel.edge);
    text = `Ребро ${sel.edge}`;
  } else if (sel?.kind === 'bar') {
    view.bars.push(sel.bar);
    text = `Стержень ${sel.bar}`;
  } else if (sel?.kind === 'finding') {
    const f = sel.finding;
    const pointToSurface = f.vertex !== undefined && f.surfaces.length > 1 && !f.points[1];
    // A vertex near another surface: that surface is the one to see.
    view.surfaces.push(...(pointToSurface ? f.surfaces.slice(-1) : f.surfaces));
    view.bars.push(...f.bars);
    if (f.vertex !== undefined) view.vertices.push(f.vertex);
    if (f.edge !== undefined) view.edges.push(f.edge);
    if (f.points.length > 1) view.segment = f.points; else view.point = f.points[0];
    text = `${KIND[f.kind] || f.kind} (${CLASS[f.class]})\n${describe(f)}${f.detail ? `\n${f.detail}` : ''}`;
  }
  viewer.highlight(view);
  $('selection').textContent = state.pending ? `${text}\n\n${state.pending.prompt}` : text;
  renderActions();
}

function button(label, action, cls = '') {
  const b = document.createElement('button');
  b.textContent = label;
  if (cls) b.className = cls;
  b.onclick = action;
  return b;
}

function renderActions() {
  const box = $('actions');
  box.innerHTML = '';
  const sel = state.selection;
  if (state.pending) {
    box.appendChild(button('Отменить выбор', () => { state.pending = null; select(sel); }));
    return;
  }
  if (!sel) return;
  if (sel.kind === 'surface') {
    box.appendChild(button('Удалить поверхность', () => edit({ op: 'delete_surface', surface: sel.surface }), 'danger'));
    box.appendChild(button('Объединить с другой…', () => waitFor('surface', 'Выберите поверхность, которая войдёт в эту', (other) =>
      edit({ op: 'join_surfaces', keep: sel.surface, other: other.surface }))));
  } else if (sel.kind === 'vertex') {
    const p = state.scene.vertices[sel.vertex];
    const row = document.createElement('div');
    row.className = 'row';
    const inputs = p.map((x) => {
      const i = document.createElement('input');
      i.type = 'number';
      i.step = '0.001';
      i.value = format(x, 4);
      row.appendChild(i);
      return i;
    });
    box.appendChild(row);
    box.appendChild(button('Сдвинуть в эти координаты', () =>
      edit({ op: 'move_vertex', vertex: sel.vertex, to: inputs.map((i) => Number(i.value)) })));
    box.appendChild(button('Слить с другой вершиной…', () => waitFor('vertex', 'Выберите вершину, которая останется', (other) =>
      edit({ op: 'merge_vertices', drop: sel.vertex, keep: other.vertex }))));
    const hint = document.createElement('div');
    hint.className = 'hint';
    hint.textContent = 'Вершина остаётся на плоскостях всех своих поверхностей; иначе правка отклоняется.';
    box.appendChild(hint);
  } else if (sel.kind === 'edge') {
    box.appendChild(button('Разбить ребро в точке клика', () => edit({ op: 'split_edge', edge: sel.edge, at: sel.point })));
  } else if (sel.kind === 'finding') {
    const f = sel.finding;
    if ((f.kind === 'gap' || f.kind === 'accepted_joint') && f.vertex !== undefined) {
      const surface = f.surfaces[f.surfaces.length - 1];
      if (f.kind === 'gap') {
        box.appendChild(button('Это шов — оставить', () => edit({ op: 'mark_joint', vertex: f.vertex, surface }, 'шов')));
      }
      box.appendChild(button('Не шов — закрыть зазор', () => edit({
        op: 'close_gap', vertex: f.vertex, surface, tolerance: f.value + 0.001,
      }, 'не шов'), 'primary'));
    }
    if (f.vertex !== undefined) {
      box.appendChild(button('Выбрать вершину', () => select({ kind: 'vertex', vertex: f.vertex })));
    }
    for (const s of f.surfaces) {
      box.appendChild(button(`Выбрать поверхность ${s}`, () => select({ kind: 'surface', surface: s })));
    }
  }
}

function waitFor(kind, prompt, done) {
  state.pending = { kind, prompt, done };
  $('pick-mode').value = kind;
  select(state.selection);
}

function onPick(event) {
  const mode = state.pending?.kind || $('pick-mode').value;
  const hit = viewer.pick(event, mode);
  if (!hit) return;
  if (state.pending) {
    const { done } = state.pending;
    state.pending = null;
    done(hit);
    return;
  }
  select({ kind: mode, ...hit });
}

async function edit(change, note = '') {
  const before = state.summary?.audit;
  const result = await busy('Правка и проверка…', async () => {
    const summary = await call('apply', { edit: change, note: note || (state.finding ? state.finding.kind : '') });
    await refresh();
    return summary;
  });
  if (result) {
    const a = result.audit;
    status(`${OP[change.op] || change.op}: ${outcome(result.last?.outcome)} · ошибок ${before?.failures ?? '?'} → ${a.failures}, PLAXIS ${before?.plaxis ?? '?'} → ${a.plaxis}`);
  }
}

/// Whether unsaved edits may be dropped: save them, discard them or stay.
async function mayDiscard() {
  if (!state.summary?.dirty) return true;
  if (await question(`Правки (${state.summary.edits}) не сохранены в проекте. Сохранить проект сейчас?`)) return saveProject();
  return question('Отбросить несохранённые правки и продолжить?');
}

async function saveProject() {
  const path = await pickFile('project', true);
  if (!path) return false;
  const saved = await busy('Сохранение проекта…', () => call('save_project', { path }));
  if (!saved) return false;
  state.summary = await call('summary');
  renderSummary();
  status(saved.input_changed_on_disk
    ? `Проект сохранён: ${path}. Внимание: файл модели на диске изменился после открытия — проект привязан к открытой версии.`
    : `Проект сохранён: ${path}`, saved.input_changed_on_disk);
  return true;
}

$('open-model').onclick = async () => {
  if (!await mayDiscard()) return;
  const path = await pickFile('model');
  if (!path) return;
  const ok = await busy('Восстановление геометрии…', async () => {
    const summary = await call('open_model', { path });
    viewer.data = null;
    await refresh();
    return summary;
  });
  if (ok) status(`Открыто: ${path}`);
};

$('open-project').onclick = async () => {
  if (!await mayDiscard()) return;
  const path = await pickFile('project');
  if (!path) return;
  const summary = await busy('Восстановление и повтор правок…', async () => {
    const s = await call('open_project', { path });
    viewer.data = null;
    await refresh();
    return s;
  });
  if (summary) {
    status(summary.replay_stopped
      ? `Повторено правок: ${summary.replayed}; остановлено: ${summary.replay_stopped}`
      : `Проект открыт, повторено правок: ${summary.replayed}${summary.input_changed ? ' (файл модели изменился)' : ''}`,
    Boolean(summary.replay_stopped));
  }
};

$('save-project').onclick = () => saveProject();

$('export-report').onclick = async () => {
  const path = await pickFile('report', true);
  if (path && await busy('Сохранение результата…', () => call('export_report', { path }))) status(`Результат сохранён: ${path}`);
};

$('undo').onclick = () => busy('Отмена…', async () => { await call('undo'); await refresh(); status('Правка отменена'); });
$('redo').onclick = () => busy('Повтор…', async () => { await call('redo'); await refresh(); status('Правка повторена'); });
$('fit').onclick = () => viewer.fit();
for (const id of ['show-failure', 'show-plaxis', 'show-review']) {
  $(id).onchange = () => { viewer.setFindings(visibleFindings()); renderFindings(); };
}
window.addEventListener('keydown', (e) => {
  if (!(e.ctrlKey || e.metaKey) || !state.summary) return;
  if (e.key === 'z' && !$('undo').disabled) { e.preventDefault(); $('undo').click(); }
  if (e.key === 'y' && !$('redo').disabled) { e.preventDefault(); $('redo').click(); }
});

// For automated checks of the interface.
window.topoEditor = { state, viewer, call, refresh, edit, select, selectFinding };
