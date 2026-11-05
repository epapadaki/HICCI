/* Draws the biological picture of one interaction as an SVG string.
   Only what the /mechanism endpoint returns is drawn; nothing is invented. */
function mechanismSvg(m, showSignal) {
  const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const W = 640, memTop = 250, memH = 34, memBot = memTop + memH;
  const H = showSignal && m.signalling.length ? 540 : 420;
  const subs = m.receptor.subunits.length ? m.receptor.subunits : [m.receptor.label];
  const co = (m.co_stimulatory || []).slice(0, 2);
  const COLORS = ["var(--accent)", "#16a34a", "#7c3aed", "#0891b2"];
  let out = [];
  const add = s => out.push(s);

  add(`<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" class="mech-svg" role="img" aria-label="${esc(m.ligand.label + " to " + m.receptor.label + ", " + m.type_label)}">`);
  add(`<defs>
    <marker id="mk-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" markerUnits="userSpaceOnUse" orient="auto"><path d="M0 0L10 5L0 10z" fill="var(--text-muted)"/></marker>
    <marker id="mk-inh" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="10" markerHeight="10" markerUnits="userSpaceOnUse" orient="auto"><path d="M5 0L5 10" stroke="var(--bad-tx)" stroke-width="2.5"/></marker>
    <pattern id="lipid" width="14" height="${memH}" patternUnits="userSpaceOnUse">
      <circle cx="7" cy="5" r="3.2" fill="var(--border-strong)"/><circle cx="7" cy="${memH - 5}" r="3.2" fill="var(--border-strong)"/>
      <path d="M5.5 8v7M8.5 8v7M5.5 ${memH - 8}v-7M8.5 ${memH - 8}v-7" stroke="var(--border)" stroke-width="1"/>
    </pattern></defs>`);

  // zones
  add(`<rect x="0" y="0" width="${W}" height="${memTop}" fill="var(--bg-soft)"/>`);
  add(`<rect x="0" y="${memBot}" width="${W}" height="${H - memBot}" fill="color-mix(in srgb, var(--accent) 4%, var(--bg-shell))"/>`);
  add(`<rect x="0" y="${memTop}" width="${W}" height="${memH}" fill="url(#lipid)"/>`);
  add(`<text x="12" y="${memTop - 8}" class="mk-zone">Extracellular</text>`);
  add(`<text x="12" y="${memBot + 18}" class="mk-zone">Receiver: ${esc(m.receiver)}</text>`);
  add(`<ellipse cx="80" cy="${H - 70}" rx="50" ry="34" fill="none" stroke="var(--border-strong)" stroke-dasharray="3 3"/><text x="80" y="${H - 66}" text-anchor="middle" class="mk-faint">nucleus</text>`);

  // receptor geometry
  const n = subs.length, gap = 112;
  const x0 = 360 - ((n - 1) * gap) / 2;
  const rx = subs.map((_, i) => x0 + i * gap);
  const headY = memTop - 58;
  const ligX = rx.reduce((a, b) => a + b, 0) / n, ligY = 120;

  // sender and ligand, by mechanism type
  const lig = m.ligand.label;
  const ligShape = (cx, cy) => m.type === "nonprotein"
    ? `<polygon points="${[0,1,2,3,4,5].map(k => `${cx + 13 * Math.cos(Math.PI / 3 * k)},${cy + 13 * Math.sin(Math.PI / 3 * k)}`).join(" ")}" fill="var(--brand-400)" stroke="var(--brand-tx)"/>`
    : `<circle cx="${cx - 9}" cy="${cy - 6}" r="11" fill="var(--brand-400)"/><circle cx="${cx + 9}" cy="${cy - 6}" r="11" fill="var(--brand-400)"/><circle cx="${cx}" cy="${cy + 9}" r="11" fill="var(--brand-500)"/>`;

  if (m.type === "contact") {
    add(`<path d="M ${ligX - 190} 40 Q ${ligX} 10 ${ligX + 190} 40 L ${ligX + 190} 120 Q ${ligX} 150 ${ligX - 190} 120 Z" fill="color-mix(in srgb, var(--brand-400) 10%, var(--bg-shell))" stroke="var(--brand-400)" stroke-width="1.5"/>`);
    add(`<text x="${ligX}" y="80" text-anchor="middle" class="mk-cell">Sender: ${esc(m.sender)}</text>`);
    add(`<rect x="${ligX - 6}" y="128" width="12" height="24" rx="3" fill="var(--brand-500)"/>`);
    add(ligShape(ligX, 162));
    add(`<text x="${ligX - 30}" y="${162}" text-anchor="end" class="mk-mol">${esc(lig)}</text><text x="${ligX - 30}" y="${182}" text-anchor="end" class="mk-faint">membrane-bound ligand</text>`);
  } else {
    add(`<ellipse cx="110" cy="100" rx="80" ry="62" fill="color-mix(in srgb, var(--brand-400) 10%, var(--bg-shell))" stroke="var(--brand-400)" stroke-width="1.5"/>`);
    add(`<ellipse cx="100" cy="104" rx="26" ry="18" fill="none" stroke="var(--brand-400)" stroke-opacity=".5"/>`);
    add(`<text x="110" y="182" text-anchor="middle" class="mk-cell">Sender: ${esc(m.sender)}</text>`);
    if (m.type === "ecm") {
      for (let k = 0; k < 4; k++) add(`<path d="M 210 ${150 + k * 16} C 300 ${120 + k * 16}, 380 ${190 + k * 16}, ${W - 30} ${150 + k * 16}" fill="none" stroke="var(--border-strong)" stroke-width="2" opacity=".6"/>`);
      add(`<text x="${W - 32}" y="140" text-anchor="end" class="mk-faint">extracellular matrix</text>`);
      add(`<path d="M 190 110 Q 240 110 270 150" fill="none" stroke="var(--text-muted)" stroke-width="1.5" marker-end="url(#mk-arr)"/><text x="232" y="100" class="mk-faint">deposited</text>`);
    } else {
      add(`<path d="M 192 104 Q ${(192 + ligX) / 2} 80 ${ligX - 34} ${ligY}" fill="none" stroke="var(--text-muted)" stroke-width="1.5" marker-end="url(#mk-arr)"/>`);
      add(`<text x="${(192 + ligX) / 2}" y="78" text-anchor="middle" class="mk-faint">${m.type === "nonprotein" ? "made and released" : "released"}</text>`);
    }
    add(ligShape(ligX, ligY));
    add(`<text x="${ligX + 26}" y="${ligY - 2}" class="mk-mol">${esc(lig)}</text><text x="${ligX + 26}" y="${ligY + 13}" class="mk-faint">${m.type === "nonprotein" ? "small molecule" : "ligand"}</text>`);
    add(`<path d="M ${ligX} ${ligY + 24} L ${ligX} ${headY - 36}" stroke="var(--text-muted)" stroke-width="1.5" marker-end="url(#mk-arr)"/>`);
  }

  // receptor subunits
  subs.forEach((g, i) => {
    const x = rx[i], col = COLORS[i % COLORS.length];
    add(`<path d="M ${x - 20} ${headY} Q ${x} ${headY + 16} ${x + 20} ${headY} L ${x + 16} ${headY + 34} L ${x - 16} ${headY + 34} Z" fill="${col}"/>`);
    add(`<rect x="${x - 7}" y="${headY + 32}" width="14" height="${memBot + 26 - headY - 32}" rx="3" fill="${col}"/>`);
    add(`<text x="${x}" y="${m.type === "contact" ? memTop - 4 : headY - 10}" text-anchor="middle" class="mk-mol" ${m.type === "contact" ? `dx="${n > 1 ? 0 : 52}"` : ""}>${esc(g)}</text>`);
  });
  if (n > 1) add(`<path d="M ${rx[0] - 26} ${memTop - 4} L ${rx[0] - 26} ${memTop + 6} L ${rx[n - 1] + 26} ${memTop + 6} L ${rx[n - 1] + 26} ${memTop - 4}" fill="none" stroke="var(--text-faint)"/><text x="${ligX}" y="${memTop + 1}" text-anchor="middle" class="mk-faint" dy="-6"></text>`);
  add(`<text x="${rx[n - 1] + 28}" y="${headY + 26}" class="mk-faint">${n > 1 ? "receptor complex" : "receptor"}</text>`);
  co.forEach((g, i) => {
    const x = rx[n - 1] + 110 + i * 60;
    add(`<path d="M ${x - 16} ${headY + 6} Q ${x} ${headY + 18} ${x + 16} ${headY + 6} L ${x + 12} ${headY + 34} L ${x - 12} ${headY + 34} Z" fill="none" stroke="#16a34a" stroke-dasharray="3 2"/>`);
    add(`<rect x="${x - 5}" y="${headY + 32}" width="10" height="${memBot + 20 - headY - 32}" rx="3" fill="none" stroke="#16a34a" stroke-dasharray="3 2"/>`);
    add(`<text x="${x}" y="${headY - 2}" text-anchor="middle" class="mk-faint">${esc(g)} (co-receptor)</text>`);
  });

  // signalling (receptor-general, SIGNOR)
  if (showSignal && m.signalling.length) {
    const top = memBot + 50;
    add(`<text x="${W - 12}" y="${H - 12}" text-anchor="end" class="mk-faint">receptor signalling in general (SIGNOR)</text>`);
    m.signalling.forEach((s, i) => {
      const fx = rx[Math.max(0, subs.indexOf(s.from))] ?? rx[0];
      const bx = 150 + i * 112, by = top + 40 + (i % 2) * 34;
      const inh = s.direction === "inhibition";
      add(`<path d="M ${fx} ${memBot + 26} Q ${fx} ${by - 26} ${bx} ${by - 17}" fill="none" stroke="${inh ? "var(--bad-tx)" : "var(--text-muted)"}" stroke-width="1.3" stroke-dasharray="4 3" marker-end="url(#${inh ? "mk-inh" : "mk-arr"})"/>`);
      add(`<rect x="${bx - 52}" y="${by - 15}" width="104" height="32" rx="${s.phenotype ? 16 : 6}" fill="var(--bg-shell)" stroke="${inh ? "var(--bad-tx)" : "var(--border-strong)"}"/>`);
      add(`<text x="${bx}" y="${by + 6}" text-anchor="middle" class="mk-step">${esc(String(s.target).replace(/_/g, " ").slice(0, 12))}</text>`);
    });
  }
  add(`</svg>`);
  return out.join("");
}
