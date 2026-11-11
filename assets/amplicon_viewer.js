(function (global) {
"use strict";

let DATA = null;
const SER = { A: 0, C: 1, G: 2, T: 3, N: 4, DEL: 5, INS: 6, SUB: 7 };
const NSER = 8;

const BASE_COL = { A: "#3da853", C: "#4285f4", G: "#f9ab00", T: "#ea4335", N: "#9aa0a6" };
const COL_DEL = "#5f6368";
const COL_REF = "#c8ccd1";
const COL_INS = "#a142f4";
const COL_SUB = "#00897b";
const COL_GUIDE = "#0e9e86";
const COL_QWIN = "#f9ab0033";

// Track layout, top to bottom. Guide rows are appended per amplicon.
const TRACKS = [
  { key: "comp", label: "Base %", h: 170 },
  { key: "mods", label: "Ins / Sub", h: 70 },
  { key: "ruler", label: "Position", h: 18 },
  { key: "ref", label: "Reference", h: 22 },
];
const GUIDE_H = 17;

// Classic scrollbars take vertical space inside the scroller and would otherwise sit
// on top of the last guide row. macOS overlay scrollbars measure 0, which is also correct.
const SBH = (() => {
  const d = document.createElement("div");
  d.style.cssText = "position:absolute;visibility:hidden;overflow-x:scroll;width:100px;height:100px";
  document.body.appendChild(d);
  const h = d.offsetHeight - d.clientHeight;
  d.remove();
  return h;
})();

let BLOB = null;          // Uint8Array, all amplicons concatenated
let amp = null;           // currently displayed amplicon record
let state = { ppb: 3, mode: "ref", hover: -1 };

const el = (id) => document.getElementById(id);
// Resolved in init(): in the report these live inside a tab that may not exist yet.
let cv, ctx, scroller, spacer, wrap;

/* ---------- data access ---------- */

// Series value at a position, as a fraction 0..1. Bytes are series-major per amplicon.
function val(series, pos) {
  return BLOB[amp.off + series * amp.len + pos] / 255;
}

function totalHeight() {
  return TRACKS.reduce((a, t) => a + t.h, 0) + amp.guides.length * GUIDE_H;
}

/* ---------- base64 + gzip ---------- */

function b64ToBytes(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

async function gunzip(bytes) {
  if (typeof DecompressionStream === "undefined") {
    throw new Error("This browser lacks DecompressionStream; use a current Chrome, Firefox, Safari or Edge.");
  }
  const ds = new DecompressionStream("gzip");
  const stream = new Blob([bytes]).stream().pipeThrough(ds);
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/* ---------- selectors ---------- */

function samplesList() {
  return [...new Set(DATA.amplicons.map((a) => a.sample))].sort();
}

function ampsFor(sample, filter) {
  const f = (filter || "").trim().toLowerCase();
  return DATA.amplicons
    .map((a, i) => ({ a, i }))
    .filter(({ a }) => a.sample === sample && (!f || a.name.toLowerCase().includes(f)));
}

function populateSamples() {
  const sel = el("av-sel-sample");
  sel.innerHTML = "";
  for (const s of samplesList()) {
    const o = document.createElement("option");
    o.value = s; o.textContent = s;
    sel.appendChild(o);
  }
}

function populateAmplicons() {
  const sel = el("av-sel-amp");
  const prev = sel.value;
  sel.innerHTML = "";
  const list = ampsFor(el("av-sel-sample").value, el("av-filter").value);
  for (const { a, i } of list) {
    const o = document.createElement("option");
    o.value = String(i);
    const flag = a.block && !a.block_ok ? " ⚠" : "";
    o.textContent = `${a.name} — ${(100 * a.modified).toFixed(1)}% mod, ${a.reads} reads${flag}`;
    sel.appendChild(o);
  }
  if (!list.length) return;
  sel.value = list.some((x) => String(x.i) === prev) ? prev : String(list[0].i);
  selectAmplicon(Number(sel.value));
}

/* ---------- header metrics ---------- */

function renderMetrics() {
  const m = el("av-metrics");
  const cutPct = (100 * amp.cut_del).toFixed(1);
  // A deletion block bounded by cut sites is the expected multi-guide outcome (excision
  // between two cuts), so it is reported, not flagged. Only a block whose ends do not
  // line up with any cut site is called out.
  const unexplained = !!(amp.block && !amp.block_ok);
  const span = amp.block
    ? `${amp.block[0]}–${amp.block[1]} (${amp.block[1] - amp.block[0] + 1} bp)`
    : "none";
  m.innerHTML = [
    ["Length", `${amp.len} bp`],
    ["Reads aligned", amp.reads.toLocaleString()],
    ["Modified", `${(100 * amp.modified).toFixed(1)}%`],
    ["Guides", String(amp.guides.length)],
    ["Deletion at cut site", `${cutPct}%`],
    [amp.block_ok ? "Excised span" : "Deletion block", span, unexplained],
  ].map(([label, value, warn]) =>
      `<div class="av-metric${warn ? " warn" : ""}"><div class="value">${value}</div><div class="label">${label}</div></div>`
    ).join("");

  const b = el("av-banner");
  if (unexplained) {
    b.className = "banner warn show";
    b.innerHTML = `<b>Deletion block ${span} does not line up with a cut site.</b>
      Inter-guide excision produces a block whose ends sit on two cut sites; this one does not,
      so the cause is worth identifying before trusting the edit call for this amplicon.`;
  } else {
    b.className = "banner warn";
  }
}

/* ---------- drawing ---------- */

function resize() {
  const w = wrap.clientWidth, h = totalHeight();
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(h * dpr);
  cv.style.height = h + "px";
  wrap.style.height = (h + SBH) + "px";
  scroller.style.height = (h + SBH) + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  draw();
}

function trackTop(key) {
  let y = 0;
  for (const t of TRACKS) {
    if (t.key === key) return y;
    y += t.h;
  }
  return y;
}

function draw() {
  if (!amp || !BLOB) return;
  const W = wrap.clientWidth, H = totalHeight(), ppb = state.ppb, sl = scroller.scrollLeft;
  ctx.clearRect(0, 0, W, H);

  const first = Math.max(0, Math.floor(sl / ppb));
  const last = Math.min(amp.len - 1, Math.ceil((sl + W) / ppb));
  const X = (p) => p * ppb - sl;

  // quantification-window shading, behind everything
  ctx.fillStyle = COL_QWIN;
  for (const p of amp.qwin) {
    if (p >= first - 1 && p <= last + 1) ctx.fillRect(X(p), 0, Math.max(1, ppb), H);
  }

  drawComposition(first, last, X, ppb);
  drawMods(first, last, X, ppb);
  drawRuler(first, last, X, ppb, W);
  drawRefSeq(first, last, X, ppb);
  drawGuides(X, ppb, H);
  drawBlock(X, ppb);
  drawCutLines(X, H);
  drawHover(X, ppb, H);

  // track separators
  ctx.strokeStyle = "#8883"; ctx.lineWidth = 1;
  let y = 0;
  for (const t of TRACKS) {
    y += t.h;
    ctx.beginPath(); ctx.moveTo(0, y + 0.5); ctx.lineTo(W, y + 0.5); ctx.stroke();
  }
}

function drawComposition(first, last, X, ppb) {
  const top = trackTop("comp"), h = TRACKS[0].h - 1;
  const w = Math.max(1, ppb - (ppb >= 4 ? 1 : 0));
  for (let p = first; p <= last; p++) {
    const x = X(p);
    const refb = amp.refseq[p];
    let acc = 0;
    const seg = (frac, col) => {
      if (frac <= 0) return;
      const hh = frac * h;
      ctx.fillStyle = col;
      ctx.fillRect(x, top + h - acc - hh, w, Math.max(0.5, hh));
      acc += hh;
    };
    if (state.mode === "ref") {
      // reference portion first (muted), then anything that is not reference
      const ri = SER[refb];
      seg(ri === undefined ? 0 : val(ri, p), COL_REF);
      for (const b of "ACGTN") {
        if (b === refb) continue;
        seg(val(SER[b], p), BASE_COL[b]);
      }
    } else {
      for (const b of "ACGTN") seg(val(SER[b], p), BASE_COL[b]);
    }
    seg(val(SER.DEL, p), COL_DEL);
  }
}

function drawMods(first, last, X, ppb) {
  const top = trackTop("mods"), h = TRACKS[1].h - 1;
  // autoscale so low-level signal stays visible, floor at 2%
  let mx = 0.02;
  for (let p = first; p <= last; p++) mx = Math.max(mx, val(SER.INS, p), val(SER.SUB, p));
  const w = Math.max(1, ppb - (ppb >= 4 ? 1 : 0));
  for (const [ser, col] of [[SER.SUB, COL_SUB], [SER.INS, COL_INS]]) {
    ctx.fillStyle = col;
    for (let p = first; p <= last; p++) {
      const v = val(ser, p) / mx;
      if (v <= 0) continue;
      const hh = Math.min(1, v) * h;
      ctx.fillRect(X(p), top + h - hh, w, Math.max(0.5, hh));
    }
  }
  ctx.fillStyle = "#8889";
  ctx.font = "10px -apple-system, Segoe UI, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`0–${(100 * mx).toFixed(mx < 0.1 ? 1 : 0)}%`, 4, top + 11);
}

function drawRuler(first, last, X, ppb, W) {
  const top = trackTop("ruler"), h = TRACKS[2].h;
  ctx.fillStyle = "#8889";
  ctx.font = "10px ui-monospace, Menlo, monospace";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  const step = ppb >= 12 ? 10 : ppb >= 6 ? 25 : ppb >= 3 ? 50 : 100;
  for (let p = Math.ceil(first / step) * step; p <= last; p += step) {
    const x = X(p) + ppb / 2;
    ctx.fillText(String(p), x, top + h / 2);
  }
}

function drawRefSeq(first, last, X, ppb) {
  const top = trackTop("ref"), h = TRACKS[3].h;
  if (ppb >= 7) {
    ctx.font = `${Math.min(13, ppb + 2)}px ui-monospace, Menlo, monospace`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    for (let p = first; p <= last; p++) {
      const b = amp.refseq[p];
      ctx.fillStyle = BASE_COL[b] || "#8889";
      ctx.fillText(b, X(p) + ppb / 2, top + h / 2);
    }
  } else {
    // too dense for letters: coloured strip
    for (let p = first; p <= last; p++) {
      ctx.fillStyle = BASE_COL[amp.refseq[p]] || "#8889";
      ctx.fillRect(X(p), top + 5, Math.max(1, ppb), h - 11);
    }
  }
}

function drawGuides(X, ppb, H) {
  const top = trackTop("ref") + TRACKS[3].h;
  amp.guides.forEach((g, gi) => {
    const y = top + gi * GUIDE_H + 3;
    const h = GUIDE_H - 7;
    const x0 = X(g.start), x1 = X(g.end + 1);
    ctx.fillStyle = COL_GUIDE;
    ctx.globalAlpha = 0.85;
    // body, with a notch at the 3' end to show strand
    const tip = Math.min(7, (x1 - x0) / 3);
    ctx.beginPath();
    if (g.strand === "+") {
      ctx.moveTo(x0, y); ctx.lineTo(x1 - tip, y); ctx.lineTo(x1, y + h / 2);
      ctx.lineTo(x1 - tip, y + h); ctx.lineTo(x0, y + h);
    } else {
      ctx.moveTo(x1, y); ctx.lineTo(x0 + tip, y); ctx.lineTo(x0, y + h / 2);
      ctx.lineTo(x0 + tip, y + h); ctx.lineTo(x1, y + h);
    }
    ctx.closePath(); ctx.fill();
    ctx.globalAlpha = 1;

    if (x1 - x0 > 68) {
      ctx.fillStyle = "#fff";
      ctx.font = "9px ui-monospace, Menlo, monospace";
      ctx.textAlign = "left"; ctx.textBaseline = "middle";
      ctx.fillText(`${g.seq} (${g.strand})`, x0 + 4, y + h / 2);
    }
  });
}

// Bracket under the base track showing the extent of the dominant deletion: green when
// its ends sit on cut sites (expected excision), amber when they do not.
function drawBlock(X, ppb) {
  if (!amp.block) return;
  const y = trackTop("mods") - 5;
  const x0 = X(amp.block[0]), x1 = X(amp.block[1] + 1);
  ctx.strokeStyle = amp.block_ok ? COL_GUIDE : "#e07a5f";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(x0, y - 4); ctx.lineTo(x0, y); ctx.lineTo(x1, y); ctx.lineTo(x1, y - 4);
  ctx.stroke();
}

function drawCutLines(X, H) {
  ctx.save();
  ctx.strokeStyle = "#ea4335cc";
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 3]);
  for (const g of amp.guides) {
    const x = Math.round(X(g.cut)) + 0.5;
    ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke();
  }
  ctx.restore();
}

function drawHover(X, ppb, H) {
  if (state.hover < 0 || state.hover >= amp.len) return;
  ctx.fillStyle = "#88888833";
  ctx.fillRect(X(state.hover), 0, Math.max(1, ppb), H);
}

/* ---------- tooltip ---------- */

function showTip(evt, pos) {
  const t = el("av-tip");
  const pct = (s) => (100 * val(s, pos)).toFixed(1);
  const refb = amp.refseq[pos];
  const rows = ["A", "C", "G", "T", "N"].map((b) => {
    const v = pct(SER[b]);
    const mark = b === refb ? " <span style='opacity:.6'>(ref)</span>" : "";
    return `<tr><td><span class="sw" style="background:${BASE_COL[b]}"></span>${b}${mark}</td><td class="n">${v}%</td></tr>`;
  }).join("");
  const guides = amp.guides
    .filter((g) => pos >= g.start && pos <= g.end)
    .map((g) => `<div style="margin-top:.25rem"><span class="sw" style="background:${COL_GUIDE}"></span>${g.seq} (${g.strand})</div>`)
    .join("");
  const atCut = amp.guides.some((g) => Math.abs(g.cut - pos) <= 1)
    ? `<div style="margin-top:.25rem;color:#ea4335">◀ cut site</div>` : "";
  const inQw = amp.qwin.includes(pos)
    ? `<div style="margin-top:.2rem;opacity:.7">in quantification window</div>` : "";

  t.innerHTML = `<b>position ${pos}</b> &middot; ref <code class="seq">${refb}</code>
    <table>${rows}
      <tr><td><span class="sw" style="background:${COL_DEL}"></span>deletion</td><td class="n">${pct(SER.DEL)}%</td></tr>
      <tr><td><span class="sw" style="background:${COL_INS}"></span>insertion</td><td class="n">${pct(SER.INS)}%</td></tr>
      <tr><td><span class="sw" style="background:${COL_SUB}"></span>substitution</td><td class="n">${pct(SER.SUB)}%</td></tr>
    </table>${guides}${atCut}${inQw}`;
  t.style.display = "block";
  const pad = 14;
  let x = evt.clientX + pad, y = evt.clientY + pad;
  const r = t.getBoundingClientRect();
  if (x + r.width > window.innerWidth - 8) x = evt.clientX - r.width - pad;
  if (y + r.height > window.innerHeight - 8) y = evt.clientY - r.height - pad;
  t.style.left = x + "px";
  t.style.top = y + "px";
}

/* ---------- gutter + legend ---------- */

function renderGutter() {
  const g = el("av-gutter");
  g.innerHTML = "";
  for (const t of TRACKS) {
    const d = document.createElement("div");
    d.style.height = t.h + "px";
    d.textContent = t.label;
    g.appendChild(d);
  }
  amp.guides.forEach((gd, i) => {
    const d = document.createElement("div");
    d.style.height = GUIDE_H + "px";
    d.textContent = `guide ${i + 1} (${gd.strand})`;
    g.appendChild(d);
  });
}

function renderLegend() {
  const items = state.mode === "ref"
    ? [[COL_REF, "reference base"], [BASE_COL.A, "A"], [BASE_COL.C, "C"], [BASE_COL.G, "G"], [BASE_COL.T, "T"], [BASE_COL.N, "N"]]
    : [[BASE_COL.A, "A"], [BASE_COL.C, "C"], [BASE_COL.G, "G"], [BASE_COL.T, "T"], [BASE_COL.N, "N"]];
  items.push([COL_DEL, "deletion"], [COL_INS, "insertion"], [COL_SUB, "substitution"], [COL_GUIDE, "guide"]);
  el("av-legend").innerHTML =
    items.map(([c, l]) => `<span><span class="sw" style="background:${c}"></span>${l}</span>`).join("") +
    `<span><span class="sw" style="background:#f9ab00"></span>quantification window</span>` +
    `<span style="opacity:.7">dashed red = cut site</span>`;
}

/* ---------- selection ---------- */

function selectAmplicon(idx) {
  amp = DATA.amplicons[idx];
  setZoom(fitZoom(), false);
  const sub = el("av-hdr-sub");
  if (sub) sub.textContent = `${DATA.meta.generated} · ${DATA.amplicons.length} amplicons across ${samplesList().length} sample(s)`;
  spacer.style.width = amp.len * state.ppb + "px";
  scroller.scrollLeft = 0;
  renderMetrics();
  renderGutter();
  renderLegend();
  resize();
}

// px/base that shows the whole amplicon, clamped to the slider's range
function fitZoom() {
  const w = wrap.clientWidth || 1200;
  return Math.max(1, Math.min(18, Math.floor((w / amp.len) * 4) / 4));
}

function setZoom(ppb, redraw) {
  state.ppb = ppb;
  el("av-zoom").value = String(ppb);
  el("av-zoom-val").textContent = Number.isInteger(ppb) ? String(ppb) : ppb.toFixed(2);
  if (redraw) {
    spacer.style.width = amp.len * ppb + "px";
    draw();
  }
}

function centreOn(pos) {
  scroller.scrollLeft = Math.max(0, pos * state.ppb - wrap.clientWidth / 2);
  draw();
}

/* ---------- events ---------- */

function bindEvents() {
el("av-sel-sample").addEventListener("change", populateAmplicons);
el("av-filter").addEventListener("input", populateAmplicons);
el("av-sel-amp").addEventListener("change", (e) => selectAmplicon(Number(e.target.value)));

el("av-zoom").addEventListener("input", (e) => {
  const centreBase = (scroller.scrollLeft + wrap.clientWidth / 2) / state.ppb;
  setZoom(Number(e.target.value), false);
  spacer.style.width = amp.len * state.ppb + "px";
  scroller.scrollLeft = Math.max(0, centreBase * state.ppb - wrap.clientWidth / 2);
  draw();
});

el("av-mode").addEventListener("click", (e) => {
  const b = e.target.closest("button[data-mode]");
  if (!b) return;
  state.mode = b.dataset.mode;
  [...el("av-mode").children].forEach((c) => c.classList.toggle("active", c === b));
  renderLegend();
  draw();
});

let cutIdx = 0;
el("av-btn-cut").addEventListener("click", () => {
  if (!amp.guides.length) return;
  cutIdx = (cutIdx + 1) % amp.guides.length;
  centreOn(amp.guides[cutIdx].cut);
});

scroller.addEventListener("scroll", draw, { passive: true });

scroller.addEventListener("mousemove", (evt) => {
  const r = scroller.getBoundingClientRect();
  const pos = Math.floor((evt.clientX - r.left + scroller.scrollLeft) / state.ppb);
  if (pos < 0 || pos >= amp.len) { state.hover = -1; el("av-tip").style.display = "none"; draw(); return; }
  state.hover = pos;
  showTip(evt, pos);
  draw();
});

scroller.addEventListener("mouseleave", () => {
  state.hover = -1;
  el("av-tip").style.display = "none";
  draw();
});

// shift+wheel is the usual horizontal-scroll gesture; make a plain wheel work too
scroller.addEventListener("wheel", (e) => {
  if (Math.abs(e.deltaY) > Math.abs(e.deltaX)) {
    scroller.scrollLeft += e.deltaY;
    e.preventDefault();
  }
}, { passive: false });

window.addEventListener("resize", () => { if (amp) resize(); });
}

/* ---------- entry points ---------- */

// Called once with the decoded payload. Safe to call when the containing tab is
// hidden: nothing is drawn until resize() sees a non-zero width.
async function init(payload) {
  DATA = payload;
  cv = el("av-cv");
  if (!cv) return;
  ctx = cv.getContext("2d");
  scroller = el("av-scroller"); spacer = el("av-spacer"); wrap = el("av-wrap");
  try {
    BLOB = await gunzip(b64ToBytes(DATA.blob));
  } catch (err) {
    cv.closest(".amplicon-viewer").innerHTML =
      `<p style="color:#c4553a">Could not decode the embedded amplicon data: ${err.message}</p>`;
    return;
  }
  bindEvents();
  populateSamples();
  populateAmplicons();
}

// A canvas inside a display:none tab measures 0 wide, so the host calls this when
// the tab becomes visible.
function onShow() {
  if (!amp || !wrap || !wrap.clientWidth) return;
  setZoom(fitZoom(), false);
  spacer.style.width = amp.len * state.ppb + "px";
  resize();
}

global.ampliconViewer = { init, onShow };

}(window));
