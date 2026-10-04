// Gap decisions in the UI: mark one gap as a joint, close another, check
// the audit counts move and the journal records both. Usage as e2e.mjs.
import { chromium } from 'playwright-core';
import path from 'node:path';

const [model, url = 'http://127.0.0.1:8787/'] = process.argv.slice(2);
const executablePath = process.env.CHROMIUM || '/opt/pw-browsers/chromium-1194/chrome-linux/chrome';
const browser = await chromium.launch({ executablePath, args: ['--use-gl=swiftshader', '--enable-unsafe-swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1400, height: 850 } });
page.on('dialog', (d) => d.accept(path.resolve(model)));
const check = (cond, what) => { if (!cond) { console.error(`FAIL: ${what}`); process.exitCode = 1; } else console.log(`ok: ${what}`); };
await page.goto(url);
await page.click('#open-model');
await page.waitForFunction(() => window.topoEditor.state.scene, null, { timeout: 30 * 60 * 1000 });
await page.waitForSelector('#busy', { state: 'hidden' });
const counts = () => page.evaluate(() => window.topoEditor.state.audit.counts);
const c0 = await counts();
check((c0.gap || 0) >= 2, `gaps found: ${c0.gap || 0}`);
// Mark the first gap as a joint.
await page.evaluate(() => {
  const f = window.topoEditor.state.audit.findings.find((x) => x.kind === 'gap');
  window.topoEditor.selectFinding(f);
});
await page.click('text=Это шов — оставить');
await page.waitForFunction(() => window.topoEditor.state.summary.edits === 1);
await page.waitForSelector('#busy', { state: 'hidden' });
const c1 = await counts();
check(c1.accepted_joint === 1 && c1.gap === c0.gap - 1, `joint accepted: gaps ${c0.gap} -> ${c1.gap}`);
// Close the next gap that closes; a refusal must name its reason and
// change nothing.
let closed = false;
let refusals = 0;
const total = c1.gap;
for (let k = 0; k < total && !closed; k++) {
  await page.evaluate((k) => {
    const f = window.topoEditor.state.audit.findings.filter((x) => x.kind === 'gap')[k];
    window.topoEditor.selectFinding(f);
  }, k);
  await page.click('text=Не шов — закрыть зазор');
  await page.waitForSelector('#busy', { state: 'hidden' });
  closed = (await page.evaluate(() => window.topoEditor.state.summary.edits)) === 2;
  if (!closed) {
    refusals++;
    const text = await page.textContent('#status');
    check(text.length > 0 && (await page.evaluate(() => document.querySelector('#status').classList.contains('error'))),
      `gap ${k} refused with a reason: ${text}`);
  }
}
const c2 = await counts();
if (closed) {
  check((c2.gap || 0) < c1.gap, `gap closed after ${refusals} refusals: ${c1.gap} -> ${c2.gap || 0}; ${await page.textContent('#status')}`);
} else {
  check((c2.gap || 0) === c1.gap, `every remaining gap refused, nothing changed (${refusals})`);
}
check((await page.evaluate(() => window.topoEditor.state.summary.audit.failures)) === 0, 'no audit failure after the edits');
const journal = await page.$$eval('#journal li', (li) => li.map((x) => x.textContent));
check(journal.length === (closed ? 2 : 1), `journal: ${journal.join(' | ')}`);
await page.screenshot({ path: '/home/user/tierc/gap.png' });
await browser.close();
