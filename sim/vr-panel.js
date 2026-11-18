// Physical widths are metres. Canvas coordinates are also the ray targets.
export const XW = 800;
export const surfaces = {
  guide: { height: 470, width: 0.50, yaw: -42, lift: -0.08 },
  controls: { height: 400, width: 0.45, yaw: 42, lift: -0.08 },
  timeline: { height: 132, width: 0.69, yaw: 0, lift: 0.40 },
};
export const TL = { x: 38, y: 79, w: 724 };
const buttons = [
  ['guide','help','Controls',535,18,229,58],
  ['controls','crater','Crater',28,79,238,112], ['controls','gulf','Gulf of Mexico',281,79,238,112], ['controls','earth','Whole Earth',534,79,238,112],
  ['controls','prev','Previous chapter',28,205,215,74], ['controls','play','Play',254,205,292,74], ['controls','next','Next chapter',557,205,215,74],
  ['controls','en','English',28,316,164,58], ['controls','hu','Magyar',202,316,164,58], ['controls','exit','Leave VR',568,316,204,58],
].map(([surface,id,label,x,y,w,h])=>({surface,id,label,x,y,w,h}));
export function panelHit(surface,x,y) {
  if (surface === 'timeline' && Math.abs(y-TL.y)<=36 && x>=TL.x-16 && x<=TL.x+TL.w+16) return 'time';
  return buttons.find(b=>b.surface===surface && x>=b.x && x<=b.x+b.w && y>=b.y && y<=b.y+b.h)?.id || 'panel';
}
const helpRows = [
  ['Trigger', 'Select or drag the timeline'], ['Left stick', 'Turn and tilt the Earth'],
  ['Grip', 'Move the world or the panel'], ['Right stick', 'Across for time, up/down for zoom'],
  ['Both grips', 'Resize and rotate the world'], ['A / X', 'Play or pause'],
  ['B / Y', 'Recenter panels'],
];
const pinchHelp = [
  ['Point and pinch', 'Select or drag the timeline'], ['Pinch the panel background', 'Move the panel'],
  ['Pinch outside the panel', 'Move the world'], ['Pinch with both hands', 'Resize and rotate the world'],
  ['Double pinch in the air', 'Recenter panels'],
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

// Three silhouettes show the change of scale before the user selects a view.
function viewGlyph(g,id,x,y,color) {
  g.save(); g.translate(x,y); g.lineWidth=2; g.strokeStyle=color;
  if(id==='earth') {
    g.rotate(-0.24);
    for(const rx of [9,21]) { g.beginPath(); g.ellipse(0,0,rx,21,0,0,Math.PI*2); g.stroke(); }
    for(const v of [-9,9]) { g.beginPath(); g.ellipse(0,v,18,4,0,0,Math.PI*2); g.stroke(); }
  } else if(id==='crater') {
    for(const [rx,ry,cy] of [[43,13,4],[33,10,1],[21,6,0],[8,3,-1]]) {
      g.beginPath(); g.ellipse(0,cy,rx,ry,-0.15,0,Math.PI*2); g.stroke();
    }
  } else {
    // Nested coast and shelf contours.
    for(let i=0;i<4;i++) {
      g.beginPath(); g.moveTo(-46,10-i*6);
      g.bezierCurveTo(-27,17-i*6,-20,-15-i*3,-3,-8-i*4);
      g.bezierCurveTo(17,0-i*4,17,21-i*5,46,2-i*4); g.stroke();
    }
  }
  g.restore();
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

export function makePanelPainter(surface, canvas, tr, chapterTicks, landmarks) {
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
  const XH = surfaces[surface].height;
  const shell = g.createLinearGradient(0,0,XW,XH);
  shell.addColorStop(0,'#1e3942'); shell.addColorStop(0.55,'#122932'); shell.addColorStop(1,'#0c1c25');
  return ({ hands, lang, preset, title, story, clock, p, playing, help, perf, hover, pressed }) => {
    g.setTransform(canvas.width/XW,0,0,canvas.height/XH,0,0); g.clearRect(0,0,XW,XH); g.textBaseline='middle';
    box(2,2,XW-4,XH-4,26,shell,'#66858b');
    if (surface==='timeline') {
      text(tr('Time'),38,35,25,sea,600); text(clock,762,35,30,ink,650,'right',625);
      box(TL.x,TL.y-3,TL.w,6,3,'#50636a');
      if(p>0) box(TL.x,TL.y-3,p*TL.w,6,3,gold);
      g.fillStyle='#607378'; for(let i=0;i<=80;i++) g.fillRect(TL.x+i*TL.w/80,TL.y+7,1,i%5===0?7:3);
      g.fillStyle=muted; for(const tick of chapterTicks) g.fillRect(TL.x+tick*TL.w,TL.y-12,2,6);
      for(const [v,label] of landmarks) text(tr(label),TL.x+v*TL.w,109,18,muted,500,v===0?'left':v===1?'right':'center');
      g.beginPath(); g.arc(TL.x+p*TL.w,TL.y,hover.includes('time')?12:9,0,2*Math.PI);
      g.fillStyle=pressed.includes('time')?ink:gold; g.fill(); return;
    }
    if(surface==='guide') {
      icon(g,'crater',53,48,gold); text('Chicxulub',84,48,27,sea,600);
      g.fillStyle='#45616a'; g.fillRect(38,92,724,1);
      g.font='500 48px Georgia, serif'; g.textAlign='left'; g.fillStyle=ink;
      g.fillText(help?tr(hands?'Hand tracking':'Controllers'):title,38,139,724);
      if(help) {
        (hands?pinchHelp:helpRows).forEach(([input,action],i)=>{
          const x=38+(i%2)*370,y=192+Math.floor(i/2)*64;
          text(tr(input),x,y,21,sea,600,'left',342); text(tr(action),x,y+26,23,ink,450,'left',342);
        });
        text(perf,760,440,18,muted,500,'right');
      } else {
        if(!storyCache.has(story)){
          let size=30,lines=storyLines(g,story,size);
          while(lines.length*(size+7)>258 && size>16) lines=storyLines(g,story,--size);
          storyCache.set(story,{size,lines});
        }
        const {size,lines}=storyCache.get(story);
        lines.forEach((line,i)=>text(line,40,195+i*(size+7),size,ink,400));
      }
    } else {
      text(tr('Views'),32,43,33,ink,650);
      text(tr(hands?'Point and pinch to select':'Point and press the trigger'),765,43,20,muted,500,'right',510);
      g.fillStyle='#45616a'; g.fillRect(30,298,740,1);
    }
    for (const b of buttons) {
      if (b.surface !== surface) continue;
      const hot = hover.includes(b.id), down = pressed.includes(b.id), primary = b.id === 'play';
      const language = b.id === 'en' || b.id === 'hu';
      const selected = b.id === lang || b.id === preset;
      const fill = selected && !down ? (hot ? '#b9eee5' : '#a1dcd3') : primary ? (down ? '#dfa456' : hot ? '#ffdaad' : gold) : down ? '#476562' : hot ? '#304f4d' : '#243638';
      box(b.x,b.y,b.w,b.h,14,fill, hot || down ? ink : primary || language || selected ? null : '#4a6060');
      const cy = b.y + b.h/2, cx = b.x + b.w/2, color = selected && !down ? '#102c2d' : primary ? '#241d12' : ink;
      const label = tr(b.id === 'play' && playing ? 'Pause' : b.id === 'help' && help ? 'Back to story' : b.label);
      const symbol = b.id === 'play' && playing ? 'pause' : b.id;
      if (['crater','gulf','earth'].includes(b.id)) {
        viewGlyph(g,b.id,cx,cy-19,selected?color:sea);
        text(label,cx,cy+25,26,color,600,'center',b.w-20);
      } else if (language || b.id === 'exit' || b.id === 'help') text(label,cx,cy,23,color,550,'center',b.w-16);
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
