// Compact low-profile VR console: a single lectern bar below eye level, an attached collapsible
// story/help card directly above it, and a tiny minimized pill for unobstructed immersion.
// Physical widths are metres; logical canvas coordinates (lw x lh) are also the ray targets.
export const surfaces = {
  bar:  { lw: 1040, lh: 276, width: 0.62, pos: [0, 0, 0],           tilt: -0.36 },
  card: { lw: 1040, lh: 264, width: 0.62, pos: [0, 0.166, -0.044],  tilt: -0.20 },
  pill: { lw: 460,  lh: 92,  width: 0.27, pos: [0, -0.045, 0.02],   tilt: -0.38 },
};
export const TL = { x: 32, y: 104, w: 976 };

const buttons = [
  // bar header utilities
  ['bar', 'story', 'Story', 684, 16, 104, 40],
  ['bar', 'help',  '?',     796, 16, 46,  40],
  ['bar', 'lang',  'EN',    850, 16, 58,  40],
  ['bar', 'min',   '',      916, 16, 46,  40],
  ['bar', 'exit',  '',      970, 16, 50,  40],
  // bar primary toolbar
  ['bar', 'prev',   '',       24,  162, 68,  64],
  ['bar', 'play',   'Play',   100, 162, 164, 64],
  ['bar', 'next',   '',       272, 162, 68,  64],
  ['bar', 'crater', 'Crater', 376, 162, 140, 64],
  ['bar', 'gulf',   'Gulf',   524, 162, 138, 64],
  ['bar', 'earth',  'Earth',  670, 162, 140, 64],
  ['bar', 'cut',    'Cut',    846, 162, 170, 64],
  // bar bottom grab handle
  ['bar', 'grab',   '',       320, 242, 400, 34],
  // card close button
  ['card', 'close', '',       972, 18,  44,  44],
  // minimized pill buttons
  ['pill', 'play',   '',      12,  12,  68,  68],
  ['pill', 'expand', 'Menu',  334, 14,  112, 64],
].map(([surface, id, label, x, y, w, h]) => ({ surface, id, label, x, y, w, h }));

export function panelHit(surface, x, y) {
  if (surface === 'bar' && Math.abs(y - TL.y) <= 32 && x >= TL.x - 16 && x <= TL.x + TL.w + 16) return 'time';
  const hit = buttons.find(b => b.surface === surface && x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h);
  if (hit) return hit.id;
  if (surface === 'bar' && y > 242) return null;
  return 'panel';
}

const helpRows = [
  ['Trigger', 'Select or drag the timeline'],        ['Left stick', 'Turn and tilt the Earth'],
  ['Drag bottom bar', 'Move the console'],           ['Right stick', 'Across for time, up/down for zoom'],
  ['Grip', 'Move the world'],                        ['A / X', 'Play or pause'],
  ['Both grips', 'Resize and rotate the world'],     ['B / Y', 'Recenter console'],
];
const pinchHelp = [
  ['Point and pinch', 'Select or drag the timeline'], ['Pinch bottom bar', 'Move the console'],
  ['Pinch outside', 'Move the world'],                ['Pinch with both hands', 'Resize and rotate the world'],
  ['Double pinch in the air', 'Recenter console'],
];

const ink = '#f4f3eb', muted = '#a9bec0', sea = '#9bded3', gold = '#ffc476';
const font = (weight, size) => weight + ' ' + size + 'px Roboto, system-ui, sans-serif';

function icon(g, id, x, y, color) {
  g.save(); g.translate(x, y); g.strokeStyle = g.fillStyle = color;
  g.lineWidth = 2.3; g.lineJoin = g.lineCap = 'round'; g.beginPath();
  switch (id) {
    case 'earth':
      g.arc(0, 0, 12, 0, 2 * Math.PI);
      g.moveTo(5, 0); g.ellipse(0, 0, 5, 12, 0, 0, 2 * Math.PI);
      g.moveTo(-11, 0); g.lineTo(11, 0);
      break;
    case 'gulf':
      for (const h of [-4.5, 4.5]) { g.moveTo(-12, h); g.bezierCurveTo(-4, h - 8, 4, h + 8, 12, h); }
      break;
    case 'crater':
      for (const [rx, ry, dy] of [[13, 6.5, 0.5], [6, 2.8, 0.5]]) { g.moveTo(rx, dy); g.ellipse(0, dy, rx, ry, -0.15, 0, 2 * Math.PI); }
      break;
    case 'cut':
      // half-cutaway block cross-section icon
      g.moveTo(-11, -5); g.lineTo(0, -10); g.lineTo(11, -5); g.lineTo(11, 5); g.lineTo(0, 10); g.lineTo(-11, 5); g.closePath();
      g.moveTo(0, -10); g.lineTo(0, 10);
      g.moveTo(0, -2); g.lineTo(11, -5);
      g.moveTo(0, 3); g.lineTo(11, 0);
      break;
    case 'prev':
      g.moveTo(-8, -8); g.lineTo(-8, 8);
      g.moveTo(8, -8); g.lineTo(-3, 0); g.lineTo(8, 8); g.closePath(); g.fill();
      break;
    case 'next':
      g.moveTo(8, -8); g.lineTo(8, 8);
      g.moveTo(-8, -8); g.lineTo(3, 0); g.lineTo(-8, 8); g.closePath(); g.fill();
      break;
    case 'play':
      g.moveTo(-6, -10); g.lineTo(10, 0); g.lineTo(-6, 10); g.closePath(); g.fill();
      break;
    case 'pause':
      g.fillRect(-7.5, -9, 5, 18); g.fillRect(2.5, -9, 5, 18);
      break;
    case 'story':
      g.moveTo(-8, -6); g.lineTo(8, -6);
      g.moveTo(-8, 0); g.lineTo(8, 0);
      g.moveTo(-8, 6); g.lineTo(3, 6);
      break;
    case 'min':
      g.moveTo(-7, -3); g.lineTo(0, 4); g.lineTo(7, -3);
      break;
    case 'expand':
      g.moveTo(-7, 3); g.lineTo(0, -4); g.lineTo(7, 3);
      break;
    case 'close':
      g.moveTo(-6, -6); g.lineTo(6, 6); g.moveTo(6, -6); g.lineTo(-6, 6);
      break;
    case 'exit':
      g.moveTo(1, -8); g.lineTo(-7, -8); g.lineTo(-7, 8); g.lineTo(1, 8);
      g.moveTo(-2, 0); g.lineTo(9, 0);
      g.moveTo(5, -4); g.lineTo(9, 0); g.lineTo(5, 4);
      break;
  }
  g.stroke(); g.restore();
}

function wrapLines(g, text, size, maxW) {
  g.font = font(400, size);
  const lines = []; let line = '';
  for (const word of text.split(/\s+/)) {
    const next = line ? line + ' ' + word : word;
    if (line && g.measureText(next).width > maxW) { lines.push(line); line = word; } else line = next;
  }
  if (line) lines.push(line);
  return lines;
}

export function makePanelPainter(surface, canvas, tr, chapterTicks, landmarks) {
  const g = canvas.getContext('2d');
  const { lw, lh } = surfaces[surface];
  const storyCache = new Map();

  const drawBox = (ctx, x, y, w, h, r, fill, stroke, sw = 2) => {
    ctx.beginPath(); ctx.roundRect(x, y, w, h, r); ctx.fillStyle = fill; ctx.fill();
    if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = sw; ctx.stroke(); }
  };
  const drawText = (ctx, s, x, y, size, color, weight = 500, align = 'left', width) => {
    ctx.font = font(weight, size); ctx.fillStyle = color; ctx.textAlign = align;
    if (width === undefined) ctx.fillText(s, x, y); else ctx.fillText(s, x, y, width);
  };

  let barBg = null, barBgKey = null;
  function ensureBarBg(state) {
    const { lang, title, preset, cut, playing, drawer, hover, pressed } = state;
    const hBtn = hover.filter(id => id && id !== 'time' && id !== 'panel').join(',');
    const pBtn = pressed.filter(id => id && id !== 'time' && id !== 'panel').join(',');
    const key = [lang, title, preset, cut, playing, drawer, hBtn, pBtn].join('|');
    if (barBg && barBgKey === key) return barBg;
    barBgKey = key;
    barBg ||= Object.assign(document.createElement('canvas'), { width: canvas.width, height: canvas.height });
    const bg = barBg.getContext('2d');
    bg.setTransform(barBg.width / lw, 0, 0, barBg.height / lh, 0, 0);
    bg.clearRect(0, 0, lw, lh); bg.textBaseline = 'middle';

    // main dark-glass shell (leaves bottom strip for the grab handle)
    const sh = bg.createLinearGradient(0, 0, lw, 244);
    sh.addColorStop(0, 'rgba(22, 46, 55, 0.96)'); sh.addColorStop(0.55, 'rgba(13, 30, 38, 0.96)'); sh.addColorStop(1, 'rgba(9, 21, 28, 0.97)');
    drawBox(bg, 2, 2, lw - 4, 238, 24, sh, 'rgba(122, 163, 171, 0.55)', 2);

    // header: gold crater glyph + brand (when Story card already shows the stage title) or live stage title
    icon(bg, 'crater', 38, 36, gold);
    drawText(bg, drawer === 'story' ? 'Chicxulub' : title, 64, 36, 24, drawer === 'story' ? sea : ink, 600, 'left', 345);

    // subtle horizontal divider under header
    bg.fillStyle = 'rgba(130, 172, 180, 0.22)'; bg.fillRect(24, 66, lw - 48, 1);

    // static timeline track, ticks, chapter notches & landmarks
    drawBox(bg, TL.x, TL.y - 3, TL.w, 6, 3, '#445a61');
    bg.fillStyle = 'rgba(155, 185, 190, 0.35)';
    for (let i = 0; i <= 80; i++) bg.fillRect(TL.x + i * TL.w / 80, TL.y + 8, 1, i % 5 === 0 ? 7 : 3);
    bg.fillStyle = muted;
    for (const tick of chapterTicks) bg.fillRect(TL.x + tick * TL.w, TL.y - 13, 2, 6);
    bg.font = font(500, 17); bg.fillStyle = muted;
    for (const [v, label] of landmarks) {
      bg.textAlign = v === 0 ? 'left' : v === 1 ? 'right' : 'center';
      bg.fillText(tr(label), TL.x + v * TL.w, 134);
    }

    // subtle vertical separators in main toolbar
    bg.fillStyle = 'rgba(130, 172, 180, 0.22)';
    bg.fillRect(358, 172, 1, 44);
    bg.fillRect(828, 172, 1, 44);

    // recessed tray behind the 3 scale presets
    drawBox(bg, 370, 156, 446, 76, 18, 'rgba(7, 17, 23, 0.55)', 'rgba(110, 148, 156, 0.25)', 1.5);

    // draw all bar buttons (except bottom grab handle)
    for (const b of buttons) {
      if (b.surface !== 'bar' || b.id === 'grab') continue;
      const hot = hover.includes(b.id), down = pressed.includes(b.id);
      const primary = b.id === 'play';
      const isExit = b.id === 'exit';
      const active = (b.id === preset) || (b.id === 'cut' && cut) || (b.id === 'story' && drawer === 'story') || (b.id === 'help' && drawer === 'help');
      const isUtil = b.y < 60;
      const r = isUtil ? 12 : 14;
      let fill = active && !down ? (hot ? '#b8efe6' : '#9edbd1')
        : primary ? (down ? '#dfa456' : hot ? '#ffdaad' : gold)
        : isExit ? (down ? '#6e3230' : hot ? '#572826' : 'rgba(92, 42, 40, 0.55)')
        : down ? '#3f5f60' : hot ? '#2c4b4e'
        : isUtil ? 'rgba(34, 56, 62, 0.72)' : '#1e3338';
      const stroke = hot || down ? ink : (primary || active) ? null : isExit ? 'rgba(235, 130, 120, 0.45)' : 'rgba(116, 152, 158, 0.35)';
      drawBox(bg, b.x, b.y, b.w, b.h, r, fill, stroke, hot || down ? 2 : 1.5);

      const cx = b.x + b.w / 2, cy = b.y + b.h / 2;
      const fg = active && !down ? '#0e282b' : primary ? '#241b0f' : isExit ? '#ffd2cc' : ink;
      const icoCol = active && !down ? '#0e282b' : primary ? '#241b0f' : isExit ? '#ff9f94' : sea;

      if (b.id === 'help') {
        drawText(bg, '?', cx, cy + 1, 23, fg, 700, 'center');
      } else if (b.id === 'lang') {
        drawText(bg, lang.toUpperCase(), cx, cy + 1, 20, fg, 700, 'center');
      } else if (b.id === 'min' || b.id === 'exit' || b.id === 'prev' || b.id === 'next') {
        icon(bg, b.id, cx, cy, icoCol);
      } else if (b.id === 'story') {
        icon(bg, 'story', b.x + 22, cy, icoCol);
        drawText(bg, tr('Story'), b.x + 38, cy + 1, 19, fg, 600, 'left', b.w - 46);
      } else {
        const label = tr(b.id === 'play' && playing ? 'Pause' : b.label);
        const sym = b.id === 'play' && playing ? 'pause' : b.id;
        const fSize = primary ? 25 : 22;
        bg.font = font(600, fSize);
        const tw = Math.min(bg.measureText(label).width, b.w - 50);
        const left = cx - (tw + 32) / 2;
        icon(bg, sym, left + 10, cy, icoCol);
        drawText(bg, label, left + 32, cy + 1, fSize, fg, 600, 'left', b.w - 50);
      }
    }

    // bottom VisionOS-style grab handle pill
    const gHot = hover.includes('grab') || pressed.includes('grab');
    drawBox(bg, 430, 255, 180, 10, 5, gHot ? ink : 'rgba(195, 224, 228, 0.45)');
    if (gHot) drawText(bg, tr('Drag to move'), 624, 260, 16, muted, 500, 'left');
    return barBg;
  }

  return state => {
    const { hands, lang, title, story, clock, p, playing, drawer, perf, hover, pressed } = state;
    g.setTransform(canvas.width / lw, 0, 0, canvas.height / lh, 0, 0);
    g.clearRect(0, 0, lw, lh); g.textBaseline = 'middle';

    if (surface === 'bar') {
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.drawImage(ensureBarBg(state), 0, 0);
      g.setTransform(canvas.width / lw, 0, 0, canvas.height / lh, 0, 0); g.textBaseline = 'middle';
      // dynamic clock & timeline progress + knob
      drawText(g, clock, 668, 36, 23, sea, 600, 'right', 250);
      if (p > 0) drawBox(g, TL.x, TL.y - 3, p * TL.w, 6, 3, gold);
      const kx = TL.x + p * TL.w, tHot = hover.includes('time'), tDown = pressed.includes('time');
      if (tHot || tDown) {
        g.beginPath(); g.arc(kx, TL.y, 18, 0, 2 * Math.PI);
        g.fillStyle = 'rgba(255, 196, 118, 0.24)'; g.fill();
      }
      g.beginPath(); g.arc(kx, TL.y, tHot || tDown ? 12 : 9, 0, 2 * Math.PI);
      g.fillStyle = tDown ? ink : gold; g.fill();
      return;
    }

    if (surface === 'card') {
      const sh = g.createLinearGradient(0, 0, lw, lh);
      sh.addColorStop(0, 'rgba(18, 39, 48, 0.94)'); sh.addColorStop(1, 'rgba(10, 23, 30, 0.95)');
      drawBox(g, 2, 2, lw - 4, lh - 4, 24, sh, 'rgba(122, 163, 171, 0.50)', 2);

      g.font = '500 32px Georgia, serif'; g.textAlign = 'left'; g.fillStyle = ink;
      g.fillText(drawer === 'help' ? tr(hands ? 'Hand tracking' : 'Controllers') : title, 36, 42, 910);

      // close button
      const cb = buttons.find(b => b.surface === 'card' && b.id === 'close');
      const cHot = hover.includes('close'), cDown = pressed.includes('close');
      drawBox(g, cb.x, cb.y, cb.w, cb.h, 12, cDown ? '#3f5f60' : cHot ? '#2c4b4e' : 'rgba(34, 56, 62, 0.65)', cHot ? ink : 'rgba(116, 152, 158, 0.35)', 1.5);
      icon(g, 'close', cb.x + cb.w / 2, cb.y + cb.h / 2, cHot ? ink : muted);

      g.fillStyle = 'rgba(130, 172, 180, 0.22)'; g.fillRect(36, 74, lw - 72, 1);

      if (drawer === 'help') {
        const rows = hands ? pinchHelp : helpRows;
        rows.forEach(([input, action], i) => {
          const col = i % 2, row = Math.floor(i / 2);
          const x = 36 + col * 494, y = 96 + row * 40;
          drawText(g, tr(input), x, y, 17, sea, 600, 'left', 460);
          drawText(g, tr(action), x, y + 19, 19, ink, 450, 'left', 460);
        });
        if (perf) drawText(g, perf, lw - 36, lh - 18, 16, muted, 500, 'right');
      } else {
        const cKey = lang + '|' + story;
        if (!storyCache.has(cKey)) {
          let size = 25, lines = wrapLines(g, story, size, lw - 72);
          while (lines.length * (size + 8) > 162 && size > 16) lines = wrapLines(g, story, --size, lw - 72);
          storyCache.set(cKey, { size, lines });
        }
        const { size, lines } = storyCache.get(cKey);
        lines.forEach((line, i) => drawText(g, line, 36, 98 + i * (size + 8), size, ink, 400));
      }
      return;
    }

    if (surface === 'pill') {
      const sh = g.createLinearGradient(0, 0, lw, lh);
      sh.addColorStop(0, 'rgba(20, 43, 52, 0.94)'); sh.addColorStop(1, 'rgba(10, 23, 30, 0.95)');
      drawBox(g, 2, 2, lw - 4, lh - 4, 44, sh, 'rgba(122, 163, 171, 0.55)', 2);

      // play/pause circular gold button
      const pb = buttons.find(b => b.surface === 'pill' && b.id === 'play');
      const pHot = hover.includes('play'), pDown = pressed.includes('play');
      drawBox(g, pb.x, pb.y, pb.w, pb.h, 34, pDown ? '#dfa456' : pHot ? '#ffdaad' : gold, pHot || pDown ? ink : null, 2);
      icon(g, playing ? 'pause' : 'play', pb.x + pb.w / 2, pb.y + pb.h / 2, '#241b0f');

      // live clock
      drawText(g, clock, 208, 46, 23, ink, 600, 'center', 228);

      // expand menu button
      const eb = buttons.find(b => b.surface === 'pill' && b.id === 'expand');
      const eHot = hover.includes('expand'), eDown = pressed.includes('expand');
      drawBox(g, eb.x, eb.y, eb.w, eb.h, 32, eDown ? '#3f5f60' : eHot ? '#2e5256' : '#1e383e', eHot ? ink : 'rgba(155, 222, 211, 0.45)', 1.5);
      icon(g, 'expand', eb.x + 26, eb.y + eb.h / 2, sea);
      drawText(g, tr('Menu'), eb.x + 44, eb.y + eb.h / 2 + 1, 20, ink, 600, 'left', eb.w - 52);
    }
  };
}
