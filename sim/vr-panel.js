// Coordinates are shared by the canvas and controller hit testing, independent of texture resolution.
export const XW = 1024, XH = 704;
export const TL = { x: 50, y: 421, w: 924 };
const buttons = [
  ['en', 'English', 374, 20, 119, 62], ['hu', 'Magyar', 497, 20, 119, 62],
  ['help', 'Controls', 638, 20, 172, 62], ['exit', 'Leave VR', 824, 20, 158, 62],
  ['prev', 'Previous chapter', 42, 471, 288, 74], ['play', 'Play', 344, 471, 336, 74], ['next', 'Next chapter', 694, 471, 288, 74],
  ['crater', 'Crater', 42, 605, 304, 74], ['gulf', 'Gulf of Mexico', 360, 605, 304, 74], ['earth', 'Whole Earth', 678, 605, 304, 74],
].map(([id, label, x, y, w, h]) => ({ id, label, x, y, w, h }));
export function panelHit(x, y) {
  if (x >= TL.x - 16 && x <= TL.x + TL.w + 16 && Math.abs(y - TL.y) <= 40) return 'time';
  return buttons.find(b => x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h)?.id || 'panel';
}

const helpRows = [
  ['Trigger', 'Select or drag the timeline'], ['Left stick', 'Turn and tilt the Earth'],
  ['Grip', 'Move the world or the panel'], ['Right stick', 'Across for time, up/down for zoom'],
  ['Both grips', 'Resize and rotate the world'], ['A / X', 'Play or pause'],
  ['B / Y', 'Bring the panel back'],
];
const ink = '#f4f3eb', muted = '#b9c9c8', sea = '#9bded3', gold = '#ffc476';
const font = (weight, size) => weight + ' ' + size + 'px Roboto, system-ui, sans-serif';

function icon(g, id, x, y, color) {
  g.save(); g.translate(x, y); g.strokeStyle = g.fillStyle = color;
  g.lineWidth = 2.3; g.lineJoin = g.lineCap = 'round'; g.beginPath();
  switch (id) {
    case 'earth': g.arc(0, 0, 13, 0, 2 * Math.PI); g.moveTo(5, 0); g.ellipse(0, 0, 5, 13, 0, 0, 2 * Math.PI); g.moveTo(-12, 1); g.lineTo(12, 1); break;
    case 'gulf': for (const h of [-5, 4]) { g.moveTo(-13, h); g.bezierCurveTo(-5, h-10, 3, h+10, 13, h); } break;
    case 'crater': for (const [rx, ry, dy] of [[14,7,0],[7,3,0]]) { g.moveTo(rx, dy); g.ellipse(0, dy, rx, ry, -0.17, 0, 2*Math.PI); } break;
    case 'prev': case 'next': { const d = id === 'next' ? 1 : -1; g.moveTo(-5*d,-8); g.lineTo(4*d,0); g.lineTo(-5*d,8); break; }
    case 'play': g.moveTo(-6,-10); g.lineTo(10,0); g.lineTo(-6,10); g.closePath(); g.fill(); break;
    case 'pause': g.fillRect(-8,-9,5,18); g.fillRect(3,-9,5,18); break;
  }
  g.stroke(); g.restore();
}

function storyLines(g, text, size) {
  g.font = font(400,size);
  const lines = []; let line = '';
  for (const word of text.split(/\s+/)) {
    const next = line ? line + ' ' + word : word;
    if (line && g.measureText(next).width > XW - 96) { lines.push(line); line = word; } else line = next;
  }
  if (line) lines.push(line);
  return lines;
}

export function makePanelPainter(canvas, tr, chapterTicks, landmarks) {
  const g = canvas.getContext('2d');
  const storyCache = new Map();
  const box = (x,y,w,h,r,fill,stroke) => {
    g.beginPath(); g.roundRect(x,y,w,h,r); g.fillStyle = fill; g.fill();
    if (stroke) { g.strokeStyle = stroke; g.lineWidth = 2; g.stroke(); }
  };
  const text = (s,x,y,size,color,weight=500,align='left',width) => {
    g.font = font(weight,size); g.fillStyle = color; g.textAlign = align;
    if (width === undefined) g.fillText(s,x,y); else g.fillText(s,x,y,width);
  };
  const shell = g.createLinearGradient(0,0,XW,XH);
  shell.addColorStop(0,'#20393e'); shell.addColorStop(0.58,'#14272d'); shell.addColorStop(1,'#102328');
  return ({ lang, preset, title, story, clock, p, playing, help, perf, hover, pressed }) => {
    g.setTransform(canvas.width / XW, 0, 0, canvas.height / XH, 0, 0); g.clearRect(0,0,XW,XH);
    g.textBaseline = 'middle';
    box(2,2,XW-4,XH-4,30,shell,'#637e81');
    box(42,29,45,45,13,'#294c50');
    icon(g,'crater',65,51,sea); text('Chicxulub',100,51,32,ink,650);
    box(369,15,252,72,20,'#112429');
    g.fillStyle = '#425455'; g.fillRect(42,94,940,1);
    text(help ? tr('Hand controls') : title,44,133,42,ink,650,'left',930);
    if (help) {
      helpRows.forEach(([input,action],i) => {
        const x = 46 + (i % 2) * 481, y = 174 + Math.floor(i/2)*52;
        text(tr(input),x,y,19,sea,650); text(tr(action),x,y+24,23,ink,450,'left',452);
      });
      text(perf,966,353,18,muted,500,'right');
    } else {
      if (!storyCache.has(story)) {
        let size = 25, lines = storyLines(g,story,size);
        while (size > 18 && lines.length * (size+5) > 205) lines = storyLines(g,story,--size);
        storyCache.set(story,{size,lines});
      }
      const {size,lines} = storyCache.get(story);
      lines.forEach((line,i) => text(line,46,181+i*(size+5),size,muted,400));
    }

    box(24,373,976,186,20,'#0e1b1e');
    text(tr('Time'),50,392,21,muted); text(clock,974,392,24,gold,600,'right',790);
    box(TL.x,TL.y-4,TL.w,8,4,'#536361');
    if (p > 0) box(TL.x,TL.y-4,TL.w*p,8,4,gold);
    g.fillStyle = '#b2bdb7'; for (const tick of chapterTicks) g.fillRect(TL.x+TL.w*tick,TL.y-9,2,5);
    for (const [position,label] of landmarks) text(tr(label),TL.x+TL.w*position,449,19,muted,500,position===0?'left':position===1?'right':'center');
    g.beginPath(); g.arc(TL.x+p*TL.w,TL.y,hover.includes('time')?13:10,0,2*Math.PI);
    g.fillStyle = pressed.includes('time') ? ink : gold; g.fill(); g.lineWidth=2; g.strokeStyle='#142020'; g.stroke();

    text(tr('Views'),44,584,24,ink,600);
    g.fillStyle = '#425b60'; g.fillRect(166,584,356,1);
    text(tr('Point and press the trigger'),978,584,20,muted,450,'right');
    for (const b of buttons) {
      const hot = hover.includes(b.id), down = pressed.includes(b.id), primary = b.id === 'play';
      const language = b.id === 'en' || b.id === 'hu';
      const selected = b.id === lang || b.id === preset;
      const fill = selected && !down ? (hot ? '#b9eee5' : '#a1dcd3') : primary ? (down ? '#dfa456' : hot ? '#ffdaad' : gold) : down ? '#476562' : hot ? '#304f4d' : '#243638';
      box(b.x,b.y,b.w,b.h,14,fill, hot || down ? ink : primary || language || selected ? null : '#4a6060');
      const cy = b.y + b.h/2, cx = b.x + b.w/2, color = selected && !down ? '#102c2d' : primary ? '#241d12' : ink;
      const label = tr(b.id === 'play' && playing ? 'Pause' : b.id === 'help' && help ? 'Back to story' : b.label);
      const symbol = b.id === 'play' && playing ? 'pause' : b.id;
      if (language || b.id === 'exit' || b.id === 'help') text(label,cx,cy,23,color,550,'center',b.w-16);
      else {
        const size = primary ? 30 : 26;
        g.font = font(600,size);
        const width = Math.min(g.measureText(label).width,b.w-66);
        const left = cx-(width+43)/2;
        icon(g,symbol,left+14,cy,primary || selected ? color : sea);
        text(label,left+43,cy,size,color,600,'left',b.w-66);
      }
    }
  };
}
