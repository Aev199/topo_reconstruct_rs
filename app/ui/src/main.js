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
  broken_bar: 'Нарушено представление стержня',
  floating_group: 'Группа висит в воздухе: нет связи с основной конструкцией',
  free_bar_end: 'Свободный конец стержня',
  lost_bar_link: 'Потеряна связь стержней: в исходной модели они были в одном узле',
  lost_surface_link: 'Потеряна связь стержня с пластиной: в исходной модели узел был общим',
  surface_not_built: 'Область КЭ не построена (нет в геометрии)',
  bar_not_built: 'Стержень не построен (нет в геометрии)',
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
  delete_bar: 'Удаление стержня',
  delete_bars: 'Удаление группы стержней',
  connect_bars: 'Общий узел стержней',
  connect_bar_to_surfaces: 'Узлы стержня с поверхностями',
  connect_surfaces: 'Общее ребро поверхностей',
};
/// Search radius of "connect" edits: the gap closure of the profile.
const CONNECT_TOLERANCE = 0.05;
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
    [/^bar piece collapsed, length ([0-9.]+)$/, (m) => `участок стержня ${m[1]} м схлопнут`],
    [/^bar deleted$/, () => 'стержень удалён'],
    [/^(\d+) bars deleted$/, (m) => `удалено стержней: ${m[1]}`],
    [/^bars share vertex (\d+), moved ([0-9.]+)$/, (m) => `общий узел ${m[1]}, конец сдвинут на ${m[2]} м`],
    [/^shared node (\d+), moved ([0-9.]+)$/, (m) => `общий узел ${m[1]}, сдвиг ${m[2]} м`],
    [/^bars share vertex (\d+)$/, (m) => `общий узел ${m[1]}`],
    [/^bar shares (\d+) node\(s\) with surfaces$/, (m) => `узлов с поверхностями: ${m[1]}`],
    [/^connected, junction ([0-9.]+)$/, (m) => `общее ребро ${m[1]} м`],
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
    ['bars_apart', 'Отклонено: стержни проходят дальше допуска друг от друга.'],
    ['skew_bars', 'Отклонено: стержни скрещиваются, не пересекаясь, — общий узел согнул бы стержень.'],
    ['parallel_bars', 'Отклонено: стержни параллельны.'],
    ['crossing_outside_bar', 'Отклонено: точка пересечения вне стержня.'],
    ['node_off_bar', 'Отклонено: узел ушёл бы со своего стержня.'],
    ['no_crossing_to_share', 'Нечего связывать: свободных пересечений не найдено.'],
    ['already_connected', 'Уже связано общим ребром.'],
    ['bar_invariant', 'Отклонено: правка нарушила бы стержень (узел вне оси, совпадающие узлы или пустой участок).'],
    ['node_beyond', 'Отклонено: узел вышел бы за соседний узел или конец стержня.'],
    ['plaxis_loader_failed', 'Загрузчик PLAXIS завершился с ошибкой — проверьте, что PLAXIS Input открыт, сервер скриптов включён, порт и пароль верны.'],
    ['no Python with plxscripting', 'Не найден Python с plxscripting — укажите путь к python.exe из поставки PLAXIS.'],
    ['no_junction', 'Поверхности не пересекаются и не примыкают.'],
    ['bar_would_collapse', 'Отклонено: стержень выродился бы.'],
    ['bar_node_on_collapsed_edge', 'Отклонено: узел стержня на схлопываемом ребре — слейте в обратном направлении.'],
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
  if (f.kind === 'floating_group') {
    const [, b, c] = f.detail.match(/bars (\d+), surfaces (\d+)/) || [];
    return `стержней ${b}, поверхностей ${c}`;
  }
  if (f.kind === 'free_bar_end') {
    return `ст. ${f.bars[0]} · верш. ${f.vertex}${f.value ? ` · ближайшее на ${format(f.value, 3)} м` : ' · рядом ничего нет'}`;
  }
  if (f.kind === 'lost_bar_link') return `ст. ${f.bars.join(' и ')} · разрыв ${format(f.value, 3)} м · ${f.detail.replace('source node', 'узел ЛИРА')}`;
  if (f.kind === 'lost_surface_link') return `ст. ${f.bars[0]} → пов. ${f.surfaces[0]} · разрыв ${format(f.value, 3)} м · ${f.detail.replace('source node', 'узел ЛИРА').replace('patch', 'область')}`;
  if (f.kind.endsWith('_not_built')) return `КЭ: ${f.value} — ${f.detail}`;
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
  $('export-plaxis').disabled = false;
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
    box.appendChild(button('Связать с другой (общее ребро)…', () => waitFor('surface', 'Выберите поверхность, пересекающую или примыкающую к этой', (other) =>
      edit({ op: 'connect_surfaces', a: sel.surface, b: other.surface }))));
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
    hint(box, 'Вершина остаётся на плоскостях всех своих поверхностей; внутренний узел стержня — на своём стержне; иначе правка отклоняется.');
  } else if (sel.kind === 'bar') {
    box.appendChild(button('Связать с другим стержнем…', () => waitFor('bar', 'Выберите стержень, пересекающий этот', (other) =>
      edit({ op: 'connect_bars', a: sel.bar, b: other.bar, tolerance: CONNECT_TOLERANCE }))));
    box.appendChild(button('Узлы в местах пересечения с поверхностями', () => edit({ op: 'connect_bar_to_surfaces', bar: sel.bar })));
    box.appendChild(button('Удалить стержень', () => edit({ op: 'delete_bar', bar: sel.bar }), 'danger'));
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
    if (f.kind === 'unshared_bar_intersection' && f.bars.length === 2) {
      box.appendChild(button('Связать стержни общим узлом', () => edit({
        op: 'connect_bars', a: f.bars[0], b: f.bars[1], tolerance: CONNECT_TOLERANCE,
      }), 'primary'));
    }
    if (f.kind === 'overlapping_bars' && f.bars.length === 2) {
      hint(box, 'Наложившиеся стержни: удалите лишний (его КЭ сохраняются в происхождении).');
      for (const b of f.bars) box.appendChild(button(`Удалить стержень ${b}`, () => edit({ op: 'delete_bar', bar: b }), 'danger'));
    }
    if ((f.kind === 'unshared_bar_surface_intersection' || f.kind === 'bar_in_surface_without_contact') && f.bars.length) {
      box.appendChild(button('Узлы стержня в пересечениях с поверхностями', () => edit({ op: 'connect_bar_to_surfaces', bar: f.bars[0] }), 'primary'));
    }
    if (f.kind.startsWith('unrepresented_') && f.surfaces.length === 2) {
      box.appendChild(button('Связать поверхности общим ребром', () => edit({
        op: 'connect_surfaces', a: f.surfaces[0], b: f.surfaces[1],
      }), 'primary'));
    }
    if (f.kind === 'coplanar_overlap' && f.surfaces.length === 2) {
      const [a, b] = f.surfaces;
      if (state.scene.surfaces[a].stiffness === state.scene.surfaces[b].stiffness) {
        box.appendChild(button('Объединить поверхности', () => edit({ op: 'join_surfaces', keep: a, other: b })));
      }
      hint(box, 'Обычно контур одной поверхности обходит вершину другой: слейте лишнюю вершину контура с вершиной соседа.');
    }
    if (f.kind === 'unshared_point_contact') {
      hint(box, 'Слейте вершину с ближайшей вершиной другой поверхности или разбейте её ребро в этой точке.');
    }
    if (f.kind === 'short_edge' || f.kind === 'short_bar_piece') {
      const ends = segmentVertices(f);
      if (ends) {
        const [a, b] = ends;
        box.appendChild(button(`Схлопнуть: оставить вершину ${a}`, () => edit({ op: 'merge_vertices', drop: b, keep: a }), 'primary'));
        box.appendChild(button(`Схлопнуть: оставить вершину ${b}`, () => edit({ op: 'merge_vertices', drop: a, keep: b })));
      }
    }
    for (const fix of f.fixes || []) {
      const d = fix.distance ? ` (${format(fix.distance, 3)} м)` : '';
      const label = {
        connect_bars: `Связать стержни ${fix.edit.a} и ${fix.edit.b}${d}`,
        merge_end_into_vertex: `Слить конец (верш. ${fix.edit.drop}) с вершиной ${fix.edit.keep}${d}`,
        move_end_onto_surface: `Посадить конец (верш. ${fix.edit.vertex}) на поверхность${d}`,
        connect_bar_to_surfaces: 'Узлы стержня с поверхностями',
        delete_group: `Удалить всю группу (${f.bars.length} ст.)`,
      }[fix.title] || fix.title;
      const danger = fix.title === 'delete_group';
      box.appendChild(button(label, () => edit(fix.edit, f.kind), danger ? 'danger' : f === sel.finding && fix === f.fixes[0] ? 'primary' : ''));
    }
    if (f.kind === 'free_bar_end' || f.kind === 'floating_group') {
      hint(box, f.kind === 'floating_group'
        ? 'Подсвечена вся группа. Если это отдельная конструкция — оставьте как есть.'
        : 'Свободный конец бывает законным (низ сваи, консоль). Вариант исправления не применяется сам — проверьте место в 3D.');
    }
    if (f.vertex !== undefined) {
      box.appendChild(button('Выбрать вершину', () => select({ kind: 'vertex', vertex: f.vertex })));
    }
    for (const b of f.bars.slice(0, 6)) {
      box.appendChild(button(`Выбрать стержень ${b}`, () => select({ kind: 'bar', bar: b })));
    }
    for (const s of f.surfaces) {
      box.appendChild(button(`Выбрать поверхность ${s}`, () => select({ kind: 'surface', surface: s })));
    }
  }
}

/// Model vertices of a finding's segment: a short edge by its id, a short
/// bar piece by its end points.
function segmentVertices(f) {
  const scene = state.scene;
  if (f.edge !== undefined) {
    const e = scene.edges.find((x) => x[0] === f.edge);
    if (e) return [e[1], e[2]];
  }
  if (f.bars.length && f.points.length > 1) {
    const at = (v, p) => Math.hypot(...scene.vertices[v].map((x, k) => x - p[k])) < 1e-6;
    const piece = scene.bars.find(([axis, a, b]) => axis === f.bars[0]
      && ((at(a, f.points[0]) && at(b, f.points[1])) || (at(b, f.points[0]) && at(a, f.points[1]))));
    if (piece) return at(piece[1], f.points[0]) ? [piece[1], piece[2]] : [piece[2], piece[1]];
  }
  return null;
}

function hint(box, text) {
  const h = document.createElement('div');
  h.className = 'hint';
  h.textContent = text;
  box.appendChild(h);
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

/// PLAXIS settings of the dialog, remembered per viewer.
const PLAXIS_FIELDS = ['plx-port', 'plx-python', 'plx-factor', 'plx-stiffness', 'plx-new', 'plx-shift'];
function loadPlaxisSettings() {
  try {
    const saved = JSON.parse(localStorage.getItem('plaxis') || '{}');
    for (const id of PLAXIS_FIELDS) {
      if (id in saved) $(id)[$(id).type === 'checkbox' ? 'checked' : 'value'] = saved[id];
    }
  } catch { /* storage unavailable */ }
}
function savePlaxisSettings() {
  try {
    localStorage.setItem('plaxis', JSON.stringify(Object.fromEntries(
      PLAXIS_FIELDS.map((id) => [id, $(id).type === 'checkbox' ? $(id).checked : $(id).value]))));
  } catch { /* storage unavailable */ }
}

$('export-plaxis').onclick = async () => {
  const a = state.summary?.audit;
  if (a && !a.passed && !await question(`Аудит геометрии не пройден (ошибок ${a.failures}). Экспортировать всё равно?`)) return;
  loadPlaxisSettings();
  const dialog = $('plaxis-dialog');
  dialog.returnValue = '';
  dialog.showModal();
  const choice = await new Promise((resolve) => dialog.addEventListener('close', () => resolve(dialog.returnValue), { once: true }));
  if (choice !== 'file' && choice !== 'run') return;
  savePlaxisSettings();
  const path = await pickFile('plaxis', true);
  if (!path) return;
  const exported = await busy('Подготовка файла обмена PLAXIS…', () => call('export_plaxis', {
    path, force_factor: Number($('plx-factor').value) || 9.80665, stiffness: $('plx-stiffness').value,
  }));
  if (!exported) return;
  const missing = exported.missing_materials.length ? `; без материала жёсткости: ${exported.missing_materials.join(', ')}` : '';
  const notes = exported.material_notes.length ? `; пересчитано материалов: ${exported.material_notes.length}` : '';
  const text = `плит ${exported.plates} (полигонов ${exported.polygons}, с отверстиями разрезано ${exported.cut_surfaces}), `
    + `балок ${exported.beams}, материалов ${exported.plate_materials}+${exported.beam_materials}${missing}${notes}`;
  if (choice === 'file') {
    status(`Файл обмена сохранён: ${path} (${text}). Загрузчик: ${exported.script}`);
    return;
  }
  const result = await busy('Построение модели в PLAXIS…', () => call('run_plaxis', {
    path,
    port: Number($('plx-port').value) || 10000,
    password: $('plx-password').value,
    python: $('plx-python').value.trim(),
    new: $('plx-new').checked,
    shift_to_origin: $('plx-shift').checked,
  }));
  if (result) {
    const r = result.report;
    const oriented = r.rectangular_beams ? `, ориентировано прямоугольных балок ${r.oriented_beams}/${r.rectangular_beams}` : '';
    status(`PLAXIS: создано плит ${r.plates ?? '?'}, балок ${r.beams ?? '?'}, материалов ${r.plate_materials ?? '?'}+${r.beam_materials ?? '?'}${oriented} за ${r.seconds ?? '?'} с (${text})`,
      r.rectangular_beams > r.oriented_beams);
  }
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
