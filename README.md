# topo_reconstruct_rs

Восстановление геотехнической геометрии (плоские поверхности и оси стержней)
из текстовой КЭ-модели ПК ЛИРА для передачи в PLAXIS и MIDAS. Цель — валидная,
связная, объяснимая геометрия, пригодная для перебивки сетки; точное
воспроизведение исходной КЭ-сетки вторично.

Состояние разработки и ближайшие задачи: [`docs/DEV_STATE.md`](docs/DEV_STATE.md).
Правила приёмки: [`docs/GEOTECHNICAL_GEOMETRY.md`](docs/GEOTECHNICAL_GEOMETRY.md).
Архитектура и требования: [`docs/RECONSTRUCTION_V2.md`](docs/RECONSTRUCTION_V2.md).
Порядок работы: [`AGENTS.md`](AGENTS.md).

## Запуск

```sh
cargo run --release -- model.txt -o report.json --mesh
```

Конвейер: распознавание осей и плоскостей → совместное решение каркаса →
топологическая сборка с геотехническими правилами → (с `--mesh`) пробная сетка.
Отчёт JSON содержит все этапы с происхождением каждой правки.

Допуски по умолчанию — профиль PLAXIS (`pipeline::Profile::plaxis()`), все
длины в единицах модели (метры):

| Параметр | По умолчанию | Смысл |
|---|---:|---|
| `--element-size` | 0.5 | целевой размер элемента |
| `--gap-closure` | 0.05 | зазоры между конструкциями закрываются (h/10) |
| `--edge-collapse` | 0.05 | короткие рёбра и стержни схлопываются |
| `--stack-offset` | 0.05 | сведение стен на ось нижней стены / в одну линию |
| `--wall-end-snap` | 0.05 | притяжка торцов стен |
| `--console-width` | 0.25 | обрезка консолей за линией стыка |
| `--crack-width` | 0.01 | трещины конвертированной сетки |
| `--simplification-cap` | 0.2 | верхняя граница допуска упрощения (½ толщины плиты) |
| `--min-opening` | 1.0 | свободные отверстия уже этого заделываются |
| `--max-opening-length` | 3.0 | отверстия длиннее этого сохраняются |
| `--keep-gap-offsets` | выкл. | не закрывать зазоры поперёк плоскости |

`--frame-cache PATH` повторно использует решённый каркас (разработка),
`TOPO_TIMING=1` печатает время этапов, `TOPO_DIAG=1` — подробности сбоев сетки.

## Проверка

```sh
cargo test
python3 scripts/check_v2_assembly.py report.json
python3 scripts/check_v2_global_geometry.py report.json --output audit.json --strict
python3 scripts/check_plaxis_profile.py report.json --output plaxis.json
python3 scripts/run_fixtures.py FIXTURE_DIR target/release/topo_reconstruct_rs OUT_DIR
```

Аудиты на Python независимы от кода реконструкции (нужны `numpy`, `shapely`:
`scripts/requirements-audit.txt`). Приватные модели и полные отчёты в
репозиторий не входят.

Прежний конвейер V1 (сшивка узлов, DXF) и режим точного воспроизведения
(`--v2-preserve-details`) удалены 2026-10-04; флаги `--v2-*` переименованы
без префикса, `--v2-mesh-preview-json PATH` → `-o PATH --mesh`.

## Редактор (Windows)

`app/` — настольное приложение: ядро на Rust в том же процессе, интерфейс на
three.js в окне WebView2 (Tauri 2). Открыть модель ЛИРА → автоматическое
восстановление → находки аудита подсвечены в 3D → правка по клику (сдвинуть /
слить вершину, удалить / объединить поверхность, разбить ребро, зазор — шов
или не шов) → каждая правка транзакционна, пишется в журнал и сразу
перепроверяется аудитом; отмена и повтор; проект хранит журнал правок;
результат сохраняется в JSON того же формата (с журналом).

```sh
cd app/ui && npm ci && npm run build            # интерфейс
cargo run --release --example editor_server     # браузерный режим: http://127.0.0.1:8787/
cd app/src-tauri && npx --prefix ../ui tauri build   # portable exe (Windows)
node app/ui/tests/e2e.mjs model.txt out/        # сквозная проверка интерфейса
```

Portable-версия (один `topo-editor.exe`, без установки) собирается вручную в
GitHub Actions (`build_app`). Кэш каркаса и данные окна хранятся в папке
`topo-editor-data` рядом с exe (если она недоступна для записи — в профиле
пользователя). Нужен WebView2 Runtime: он входит в Windows 11 и обновлённую
Windows 10.
