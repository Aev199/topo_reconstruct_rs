// End-to-end check of the editor UI against the real Rust core through the
// development bridge. Start `cargo run --release --example editor_server`
// first. Usage: node tests/e2e.mjs MODEL.txt OUT_DIR [URL]
import { chromium } from 'playwright-core';
import fs from 'node:fs';
import path from 'node:path';

const [model, out, url = 'http://127.0.0.1:8787/'] = process.argv.slice(2);
if (!model || !out) {
  console.error('usage: node tests/e2e.mjs MODEL.txt OUT_DIR [URL]');
  process.exit(2);
}
fs.mkdirSync(out, { recursive: true });
const executablePath = process.env.CHROMIUM || '/opt/pw-browsers/chromium-1194/chrome-linux/chrome';
const browser = await chromium.launch({ executablePath, args: ['--use-gl=swiftshader', '--enable-unsafe-swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1400, height: 850 } });
const errors = [];
page.on('pageerror', (e) => errors.push(String(e)));
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
const answers = [];
const asked = [];
page.on('dialog', (d) => {
  const answer = answers.shift();
  if (d.type() === 'confirm') asked.push(d.message());
  return answer === false ? d.dismiss() : d.accept(answer ?? '');
});
const check = (cond, what) => {
  if (!cond) { console.error(`FAIL: ${what}`); process.exitCode = 1; } else console.log(`ok: ${what}`);
};
const summary = () => page.evaluate(() => window.topoEditor.state.summary);

await page.goto(url);
answers.push(path.resolve(model));
await page.click('#open-model');
await page.waitForFunction(() => window.topoEditor.state.scene, null, { timeout: 30 * 60 * 1000 });
await page.waitForSelector('#busy', { state: 'hidden' });
let s = await summary();
check(s.surfaces > 0, `model opened: ${s.surfaces} surfaces, ${s.bars} bars`);
const header = await page.textContent('#audit-summary');
check(header.includes('Геометрия и связность') && header.includes('Профиль PLAXIS'), `two verdicts: ${header}`);
const findings = await page.$$eval('#findings li', (li) => li.length);
check(findings === (await page.evaluate(() => window.topoEditor.state.audit.findings.filter((f) => f.class !== 'review').length)),
  `findings listed: ${findings}`);
await page.screenshot({ path: path.join(out, '1-opened.png') });
// The canvas is drawn (not a blank clear colour).
const drawn = await page.evaluate(() => {
  const c = document.querySelector('#view canvas');
  const g = c.getContext('webgl2') || c.getContext('webgl');
  const px = new Uint8Array(4 * 100);
  g.readPixels(Math.floor(c.width / 2) - 50, Math.floor(c.height / 2), 100, 1, g.RGBA, g.UNSIGNED_BYTE, px);
  return new Set(Array.from({ length: 100 }, (_, i) => px.slice(i * 4, i * 4 + 3).join())).size;
});
check(drawn > 1, `geometry drawn (${drawn} colours across the centre line)`);

if (findings) {
  await page.click('#findings li');
  const sel = await page.textContent('#selection');
  check(sel.length > 0 && sel !== 'ничего', `finding selected: ${sel.split('\n')[0]}`);
  await page.screenshot({ path: path.join(out, '2-finding.png') });
}

// Pick a surface by clicking the centre of the view in surface mode.
await page.selectOption('#pick-mode', 'surface');
const box = await page.locator('#view canvas').boundingBox();
await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
let picked = await page.evaluate(() => window.topoEditor.state.selection);
if (!picked || picked.kind !== 'surface') {
  await page.evaluate(() => window.topoEditor.select({ kind: 'surface', surface: 0 }));
  picked = await page.evaluate(() => window.topoEditor.state.selection);
  console.log('note: centre click hit no surface; selected surface 0 directly');
}
check(picked?.kind === 'surface', `surface ${picked?.surface} selected`);
const before = await summary();
await page.click('text=Удалить поверхность');
await page.waitForFunction((n) => window.topoEditor.state.summary.edits === n + 1, before.edits);
await page.waitForSelector('#busy', { state: 'hidden' });
s = await summary();
check(s.surfaces === before.surfaces - 1, `surface deleted: ${before.surfaces} -> ${s.surfaces}, audit failures ${before.audit.failures} -> ${s.audit.failures}`);
await page.screenshot({ path: path.join(out, '3-deleted.png') });
await page.click('#undo');
await page.waitForFunction((n) => window.topoEditor.state.summary.edits === n, before.edits);
await page.waitForSelector('#busy', { state: 'hidden' });
s = await summary();
check(s.surfaces === before.surfaces && s.can_redo, 'undo restores the surface');
await page.click('#redo');
await page.waitForFunction((n) => window.topoEditor.state.summary.edits === n + 1, before.edits);
await page.waitForSelector('#busy', { state: 'hidden' });
check((await summary()).surfaces === before.surfaces - 1, 'redo deletes it again');

// Project save / reopen replays the journal; the edited result is saved.
const project = path.resolve(out, 'edited.topo.json');
answers.push(project);
await page.click('#save-project');
await page.waitForFunction(() => document.querySelector('#status').textContent.includes('Проект сохранён'));
check(fs.existsSync(project), 'project saved');
answers.push(project);
await page.click('#open-project');
await page.waitForFunction(() => document.querySelector('#status').textContent.includes('Проект открыт')
  || document.querySelector('#status').classList.contains('error'), null, { timeout: 30 * 60 * 1000 });
const statusText = await page.textContent('#status');
check(statusText.includes('повторено правок: 1'), `project reopened: ${statusText}`);
check((await summary()).dirty === false, 'reopened project has no unsaved edits');
const journalText = await page.textContent('#journal');
check(journalText.includes('Удаление поверхности'), `journal in Russian: ${journalText}`);
// An unsaved edit: opening another project asks first; "no, no" cancels.
await page.click('#undo');
await page.waitForFunction(() => window.topoEditor.state.summary.dirty === true);
await page.waitForSelector('#busy', { state: 'hidden' });
asked.length = 0;
answers.push(false, false);
await page.click('#open-project');
await page.waitForFunction(() => true);
await new Promise((r) => setTimeout(r, 500));
s = await summary();
check(asked.length === 2 && s.dirty === true && s.edits === 0, `unsaved edits kept after cancel (${asked.length} questions)`);
await page.click('#redo');
await page.waitForFunction(() => window.topoEditor.state.summary.dirty === false);
await page.waitForSelector('#busy', { state: 'hidden' });
const report = path.resolve(out, 'edited-report.json');
answers.push(report);
await page.click('#export-report');
await page.waitForFunction(() => document.querySelector('#status').textContent.includes('Результат сохранён'));
check(fs.existsSync(report), 'edited result saved');
check(errors.length === 0, `no page errors${errors.length ? `: ${errors.join(' | ')}` : ''}`);
await browser.close();
