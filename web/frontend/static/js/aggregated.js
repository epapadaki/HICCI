/* Aggregated results: the curated network, built only from interactions kept
   after expert review. Same circular SVG layout and matrix as the Network page. */
function aggregatedApp() {
  const NS = "http://www.w3.org/2000/svg";
  return Object.assign(docPanel(), {
    data: null, view: "network", maxEdges: 80, typesOpen: false,
    section: "network", impact: null, _ireq: 0, hubs: null, _hreq: 0,
    profiles: null, _preq: 0,
    bio: null, _breq: 0, bioSource: "reactome", bioSel: null, bioDetail: null, _bdreq: 0,
    evid: null, _ereq: 0, evLoading: false, evView: "upset",
    enrich: null, _fereq: 0, feLoading: false,
    drill: null, _dreq: 0,
    insights: null, _nreq: 0,
    modalOpen: false,
    cutoff: 95,
    f: { types: [] },
    sel: null, pairData: null, _req: 0, _counts: {},

    async init() {
      try { const c = parseFloat(localStorage.getItem("ccc-keep-cutoff")); if (!isNaN(c)) this.cutoff = c; } catch (e) {}
      this.initPanel();
      await this.load();
      window.addEventListener("resize", () => { this.draw(); this.drawSankey(); this.drawDotplot(); });
    },

    params() {
      const p = new URLSearchParams({ cutoff: this.cutoff });
      this.f.types.forEach(t => p.append("cell_type", t));
      return p;
    },
    setCutoff(v) {
      let c = parseFloat(v); if (isNaN(c)) c = 95; c = Math.max(0, Math.min(100, c));
      this.cutoff = c;
      try { localStorage.setItem("ccc-keep-cutoff", String(c)); } catch (e) {}
      this.closeDrill();
      this.load();
      if (this.impact !== null) this.loadImpact();
      if (this.hubs !== null) this.loadHubs();
      if (this.profiles !== null) this.loadProfiles();
      if (this.bio !== null) this.loadBio();
      if (this.evid !== null) this.loadEvidence();
      if (this.enrich !== null) this.loadEnrich();
    },

    // ---------------------------------------------------------- HITL impact
    // ---- Ligand-receptor rewiring (replaces the old HITL-impact networks) ----
    async loadImpact() {
      const req = ++this._ireq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/lr-rewiring?cutoff=" + this.cutoff);
        if (req !== this._ireq) return;
        this.impact = d;
        this.$nextTick(() => this.drawRewiring());
      } catch (e) {}
    },
    _famColor(pw) {
      // Muted, harmonious categorical palette (desaturated — reads as a figure, not a toy).
      const FAM = { MIF: "#4E7C8A", PVR: "#9C6B72", CD86: "#9C6B72", GALECTIN: "#9C6B72",
        CXCL: "#5E7E9E", CCL: "#5E7E9E", TGFb: "#B08D57", VEGF: "#9A6C86",
        COMPLEMENT: "#7A9E7E", CLEC: "#7B6D8D", SPP1: "#78896B", CD45: "#5E8C9E" };
      if (FAM[pw]) return FAM[pw];
      const PAL = ["#4E7C8A", "#C08552", "#7A9E7E", "#9C6B72", "#7B6D8D", "#5E8C9E",
        "#B08D57", "#8E7A94", "#6A8E8A", "#A5785E", "#78896B", "#8C6E82"];
      let h = 0; for (const c of String(pw)) h = (h * 31 + c.charCodeAt(0)) >>> 0;
      return PAL[h % PAL.length];
    },
    // Generic multi-stage flow into one <svg>. stages = [{label, items:[names]}];
    // links = [{si, s, t, w, color}]. Node height = max(incoming, outgoing) sum.
    _flow(svgId, stages, links) {
      const svg = document.getElementById(svgId);
      if (!svg) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const NS2 = NS, el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS2, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };
      const FONT = "Inter, 'Segoe UI', system-ui, sans-serif";
      const BASE_OP = 0.4;
      const nStage = stages.length;
      // per-node incoming/outgoing sums
      const outS = stages.map(() => ({})), inS = stages.map(() => ({}));
      links.forEach(lk => { outS[lk.si][lk.s] = (outS[lk.si][lk.s] || 0) + lk.w; inS[lk.si + 1][lk.t] = (inS[lk.si + 1][lk.t] || 0) + lk.w; });
      const tot = stages.map((st, si) => st.items.map(n => Math.max(outS[si][n] || 0, inS[si][n] || 0)));
      const maxItems = Math.max(...stages.map(s => s.items.length), 1);
      const rowH = 32, pad = 40, nodeW = 12, gap = 10;
      const H = pad * 2 + 28 + rowH * maxItems;
      const W = 960, leftM = 172, rightM = 224;
      svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
      el("rect", { x: 0, y: 0, width: W, height: H, fill: "#ffffff" });
      const xs = stages.map((_, i) => leftM + (nStage > 1 ? i * ((W - leftM - rightM) / (nStage - 1)) : 0));
      const avail = H - pad * 2 - 24;
      // One shared px-per-interaction across columns (so flows stay conserved),
      // capped so a few interactions don't become crayon-thick bands, and each
      // column centred vertically to leave clean whitespace.
      const sums = tot.map(a => a.reduce((x, y) => x + y, 0) || 1);
      const maxSum = Math.max(...sums, 1);
      let unit = (avail - gap * (maxItems - 1)) / maxSum;
      unit = Math.min(unit, 13);
      const pos = {};
      stages.forEach((st, si) => {
        const colH = tot[si].reduce((a, b) => a + Math.max(12, b * unit), 0) + gap * (st.items.length - 1);
        let y = pad + 24 + Math.max(0, (avail - colH) / 2);
        pos[si] = {};
        st.items.forEach((n, k) => { const h = Math.max(12, tot[si][k] * unit); pos[si][n] = { y, h, off: 0, offR: 0, oS: outS[si][n] || 1, iS: inS[si][n] || 1, val: tot[si][k] }; y += h + gap; });
      });
      // stage headers (uppercase, tracked) with a thin rule under each column
      stages.forEach((st, si) => {
        el("text", { x: xs[si] + nodeW / 2, y: 20, "text-anchor": "middle", "font-family": FONT,
          "font-size": 11, "font-weight": 700, "letter-spacing": "0.08em", fill: "#8a8f98" }).textContent = st.label.toUpperCase();
      });
      // ribbons
      const ribbons = [];
      links.forEach(lk => {
        const a = pos[lk.si][lk.s], b = pos[lk.si + 1][lk.t];
        if (!a || !b) return;
        const uA = a.h / a.oS, uB = b.h / b.iS;
        const y1 = a.y + a.offR, y1b = y1 + lk.w * uA; a.offR += lk.w * uA;
        const y2 = b.y + b.off, y2b = y2 + lk.w * uB; b.off += lk.w * uB;
        const x1 = xs[lk.si] + nodeW, x2 = xs[lk.si + 1], mx = (x1 + x2) / 2;
        const path = el("path", { d: `M${x1},${y1.toFixed(1)} C${mx},${y1.toFixed(1)} ${mx},${y2.toFixed(1)} ${x2},${y2.toFixed(1)} L${x2},${y2b.toFixed(1)} C${mx},${y2b.toFixed(1)} ${mx},${y1b.toFixed(1)} ${x1},${y1b.toFixed(1)} Z`,
          fill: lk.color, opacity: BASE_OP, stroke: lk.color, "stroke-width": 0.5, "stroke-opacity": 0.4 });
        const title = el("title", {}, path); title.textContent = `${lk.s} → ${lk.t}`;
        ribbons.push({ path, si: lk.si, s: lk.s, t: lk.t });
      });
      // nodes + labels
      const nodeFill = ["#C07A53", "#4E7C8A", "#8A8F98"];
      const trunc = (s) => { s = String(s); return s.length > 30 ? s.slice(0, 29) + "…" : s; };
      const nodeEls = [];
      stages.forEach((st, si) => {
        const last = si === nStage - 1;
        const fill = nodeFill[si] || "#64748b";
        st.items.forEach(n => {
          const p = pos[si][n];
          const rect = el("rect", { x: xs[si], y: p.y.toFixed(1), width: nodeW, height: p.h.toFixed(1), rx: 3,
            fill, stroke: "#ffffff", "stroke-width": 1 });
          const rt = el("title", {}, rect); rt.textContent = `${n} · ${p.val}`;
          const lx = si === 0 ? xs[si] - 8 : xs[si] + nodeW + 8;
          const anchor = si === 0 ? "end" : "start";
          const label = el("text", { x: lx, y: (p.y + p.h / 2 + 4).toFixed(1), "text-anchor": anchor, "font-family": FONT, "font-size": 12, "font-weight": 600,
            fill: "#2a2d34", stroke: (si > 0 && !last) ? "#ffffff" : "none", "stroke-width": 3, "paint-order": "stroke" });
          label.textContent = trunc(n);
          const cnt = el("tspan", { fill: "#9aa0ad", "font-weight": 500, "font-size": 10.5 }, label);
          cnt.textContent = "  " + p.val;
          nodeEls.push({ si, n, rect, label });
        });
      });
      // hover: emphasise ribbons touching a node, dim the rest
      const setFocus = (si, n) => {
        ribbons.forEach(r => {
          const hit = n == null || (r.si === si && r.s === n) || (r.si + 1 === si && r.t === n);
          r.path.setAttribute("opacity", n == null ? BASE_OP : (hit ? 0.9 : 0.07));
        });
        nodeEls.forEach(ne => { ne.label.setAttribute("opacity", n == null || ne.n === n ? 1 : 0.35); });
      };
      nodeEls.forEach(ne => {
        ne.rect.style.cursor = "pointer";
        ne.rect.addEventListener("mouseenter", () => setFocus(ne.si, ne.n));
        ne.rect.addEventListener("mouseleave", () => setFocus(null, null));
      });
    },
    drawRewiring() {
      const d = this.impact; if (!d || !d.pairs) return;
      const pairs = d.pairs;
      const sumBy = (keyFn, rows) => { const m = {}; (rows || pairs).forEach(p => { const k = keyFn(p); m[k] = (m[k] || 0) + p.count; }); return m; };
      const ord = (m) => Object.keys(m).sort((a, b) => m[b] - m[a]);
      const ligs = ord(sumBy(p => p.ligand)), rcps = ord(sumBy(p => p.receptor));
      // Method-independent colour: the receptor's top Reactome pathway (falls back to
      // the receptor name so every link still gets a stable colour).
      const rl = d.reactome_links || [];
      const recTop = {};
      rl.forEach(x => { if (!recTop[x.receptor] || x.weight > recTop[x.receptor].w) recTop[x.receptor] = { pw: x.pathway, w: x.weight }; });
      const lrColor = p => this._famColor((recTop[p.receptor] && recTop[p.receptor].pw) || p.receptor);
      const lrLinks = pairs.map(p => ({ si: 0, s: p.ligand, t: p.receptor, w: p.count, color: lrColor(p) }));
      // 1) ligand -> receptor (no pathway)
      this._flow("lr-bip", [{ label: "Ligand", items: ligs }, { label: "Receptor", items: rcps }], lrLinks);
      // 2) ligand -> receptor -> Reactome pathway
      if (rl.length) {
        const trunc = (s) => s.length > 32 ? s.slice(0, 31) + "…" : s;
        const rec = new Set(rl.map(x => x.receptor));
        const rp = pairs.filter(p => rec.has(p.receptor));
        const reLigs = ord(sumBy(p => p.ligand, rp)), reRcps = ord(sumBy(p => p.receptor, rp));
        const rePws = (d.reactome_pathways || []).map(trunc);
        const reLinks = rp.map(p => ({ si: 0, s: p.ligand, t: p.receptor, w: p.count, color: lrColor(p) }))
          .concat(rl.map(x => ({ si: 1, s: x.receptor, t: trunc(x.pathway), w: x.weight, color: this._famColor(x.pathway) })));
        this._flow("lr-tri-re", [{ label: "Ligand", items: reLigs }, { label: "Receptor", items: reRcps }, { label: "Reactome pathway", items: rePws }], reLinks);
      }
    },
    async loadHubs() {
      const req = ++this._hreq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/molecular-hubs?cutoff=" + this.cutoff);
        if (req !== this._hreq) return;
        this.hubs = d;
      } catch (e) {}
    },
    hubBar(v, max) { return Math.max(3, Math.round((v / (max || 1)) * 100)) + "%"; },
    async loadProfiles() {
      const req = ++this._preq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/communication-profiles?cutoff=" + this.cutoff);
        if (req !== this._preq) return;
        this.profiles = d;
      } catch (e) {}
    },
    totalBar(v, max) { return Math.max(4, Math.round((v / (max || 1)) * 100)) + "%"; },
    sendPct(p) { return Math.round((p.balance || 0) * 100); },

    // -------------------------------------------------- biological processes
    async loadBio() {
      const req = ++this._breq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/biological-processes?cutoff=" + this.cutoff);
        if (req !== this._breq) return;
        this.bio = d; this.bioSel = null; this.bioDetail = null;
        this.$nextTick(() => this.drawSankey());
      } catch (e) {}
    },
    bioCur() { return this.bio ? this.bio[this.bioSource] : null; },
    setBioSource(s) { this.bioSource = s; this.bioSel = null; this.bioDetail = null; this.$nextTick(() => this.drawSankey()); },
    procBar(v, max) { return Math.max(3, Math.round((v / (max || 1)) * 100)) + "%"; },
    pctLabeled() { const c = this.bioCur(); return c && c.total_kept ? Math.round((c.labeled / c.total_kept) * 100) : 0; },
    async selectProcess(name) {
      this.bioSel = name; this.bioDetail = null;
      const req = ++this._bdreq;
      try {
        const p = new URLSearchParams({ source: this.bioSource, name, cutoff: this.cutoff });
        const d = await getJSON("/api/biological-processes/detail?" + p);
        if (req !== this._bdreq) return;
        this.bioDetail = d;
      } catch (e) {}
    },
    _trunc(s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + "…" : s; },
    drawSankey() {
      const cur = this.bioCur();
      const svg = document.getElementById("bio-sankey");
      if (!svg || !cur) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const sk = cur.sankey || {};
      const P = (sk.processes || []).filter(d => d.value > 0).slice(0, 10);
      if (!P.length) return;
      const pset = new Set(P.map(d => d.name));
      const linksSP = (sk.links_sp || []).filter(l => pset.has(l.target));
      const linksPR = (sk.links_pr || []).filter(l => pset.has(l.source));
      const sVal = {}, rVal = {};
      linksSP.forEach(l => sVal[l.source] = (sVal[l.source] || 0) + l.value);
      linksPR.forEach(l => rVal[l.target] = (rVal[l.target] || 0) + l.value);
      const S = Object.keys(sVal).map(n => ({ name: n, value: sVal[n] })).sort((a, b) => b.value - a.value);
      const R = Object.keys(rVal).map(n => ({ name: n, value: rVal[n] })).sort((a, b) => b.value - a.value);
      const total = P.reduce((a, d) => a + d.value, 0) || 1;

      const W = 820, maxN = Math.max(S.length, P.length, R.length, 1);
      const H = Math.max(340, 40 * maxN + 40);
      svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
      const el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };

      const colX = { s: 140, p: 404, r: 668 }, nodeW = 13, pad = 30, gap = 9;
      const scale = (H - 2 * pad - gap * (maxN - 1)) / total;
      const layout = (list, x) => { let y = pad; const pos = {}; list.forEach(d => { const h = Math.max(2, d.value * scale); pos[d.name] = { x, y0: y, h }; y += h + gap; }); return pos; };
      const posS = layout(S, colX.s), posP = layout(P, colX.p), posR = layout(R, colX.r);
      const sIdx = {}, pIdx = {}, rIdx = {};
      S.forEach((d, i) => sIdx[d.name] = i); P.forEach((d, i) => pIdx[d.name] = i); R.forEach((d, i) => rIdx[d.name] = i);

      const gL = el("g", {});
      const ribbon = (x0, ya, x1, yb, w) => { const xm = (x0 + x1) / 2; return `M ${x0} ${ya} C ${xm} ${ya} ${xm} ${yb} ${x1} ${yb} L ${x1} ${yb + w} C ${xm} ${yb + w} ${xm} ${ya + w} ${x0} ${ya + w} Z`; };
      const drawSide = (links, xa, xb, posA, idxA, posB, idxB, color) => {
        const offA = {}, offB = {};
        Object.keys(posA).forEach(k => offA[k] = posA[k].y0);
        Object.keys(posB).forEach(k => offB[k] = posB[k].y0);
        const bySrc = links.slice().sort((a, b) => idxA[a.source] - idxA[b.source] || idxB[a.target] - idxB[b.target]);
        const srcY = new Map(); bySrc.forEach(l => { srcY.set(l, offA[l.source]); offA[l.source] += l.value * scale; });
        const byTgt = links.slice().sort((a, b) => idxB[a.target] - idxB[b.target] || idxA[a.source] - idxA[b.source]);
        const tgtY = new Map(); byTgt.forEach(l => { tgtY.set(l, offB[l.target]); offB[l.target] += l.value * scale; });
        links.forEach(l => {
          const w = Math.max(1, l.value * scale);
          const p = el("path", { d: ribbon(xa, srcY.get(l), xb, tgtY.get(l), w), fill: color, "fill-opacity": 0.22, class: "sankey-link" }, gL);
          const tt = el("title", {}, p); tt.textContent = `${l.source} → ${l.target}: ${l.value}`;
        });
      };
      drawSide(linksSP, colX.s + nodeW, colX.p, posS, sIdx, posP, pIdx, "var(--accent)");
      drawSide(linksPR, colX.p + nodeW, colX.r, posP, pIdx, posR, rIdx, "var(--layer-data)");

      const gN = el("g", {});
      const drawNodes = (list, pos, fill, labelX, anchor, cls, onClick) => {
        list.forEach(d => {
          const p = pos[d.name];
          const g = el("g", { class: onClick ? "sankey-node" : "" }, gN);
          el("rect", { x: p.x, y: p.y0, width: nodeW, height: p.h, rx: 2, fill }, g);
          const t = el("text", { x: labelX(p), y: p.y0 + p.h / 2 + 3.5, "text-anchor": anchor, class: "sankey-label " + cls }, g);
          t.textContent = this._trunc(d.name, 24);
          const tt = el("title", {}, g); tt.textContent = `${d.name}: ${d.value}`;
          if (onClick) { g.addEventListener("click", () => onClick(d.name)); el("rect", { x: p.x - 2, y: p.y0 - 2, width: nodeW + 4, height: p.h + 4, fill: "transparent", style: "cursor:pointer" }, g); }
        });
      };
      drawNodes(S, posS, "var(--accent)", p => p.x - 6, "end", "", null);
      drawNodes(P, posP, "var(--text-muted)", p => p.x + nodeW + 6, "start", "p", n => this.selectProcess(n));
      drawNodes(R, posR, "var(--layer-data)", p => p.x + nodeW + 6, "start", "", null);

      for (const [x, lbl, anch] of [[colX.s + nodeW / 2, "Sender", "middle"], [colX.p + nodeW / 2, "Process", "middle"], [colX.r + nodeW / 2, "Receiver", "middle"]])
        el("text", { x, y: 16, "text-anchor": anch, class: "sankey-head" }).textContent = lbl;
    },

    // ----------------------------------------------- drill-down to CCIs
    async hubDrill(kind, name) {
      const title = (kind === "ligand" ? "Ligand " : "Receptor ") + name;
      const req = ++this._dreq; this.drill = { scope: "hubs", mode: "list", title, rows: [], pairs: [], total: 0 };
      try { const p = new URLSearchParams({ kind, name, cutoff: this.cutoff });
        const d = await getJSON("/api/molecular-hubs/detail?" + p); if (req !== this._dreq) return;
        this.drill = { scope: "hubs", mode: "list", title, ...d };
      } catch (e) {}
    },
    async hitlDrill(bucket) {
      const names = { accepted: "Explicitly accepted", rejected: "Explicitly rejected", unreviewed: "Unreviewed", kept: "Kept set" };
      const title = names[bucket] || bucket;
      const req = ++this._dreq; this.drill = { scope: "impact", mode: "list", title, rows: [], pairs: [], total: 0 };
      try { const p = new URLSearchParams({ bucket, cutoff: this.cutoff });
        const d = await getJSON("/api/hitl-impact/detail?" + p); if (req !== this._dreq) return;
        this.drill = { scope: "impact", mode: "list", title, ...d };
      } catch (e) {}
    },
    async profileDrill(population) {
      const req = ++this._dreq; this.drill = { scope: "profiles", mode: "profile", title: population, outgoing: { total: 0, rows: [] }, incoming: { total: 0, rows: [] } };
      try { const p = new URLSearchParams({ population, cutoff: this.cutoff });
        const d = await getJSON("/api/communication-profiles/detail?" + p); if (req !== this._dreq) return;
        this.drill = { scope: "profiles", mode: "profile", title: population, ...d };
      } catch (e) {}
    },
    async evidenceDrill(source) {
      const s = this.evid && this.evid.sources ? this.evid.sources.find(x => x.key === source) : null;
      if (s && !s.supported) return;
      const title = (s ? s.label : source) + " — supported interactions";
      const req = ++this._dreq; this.drill = { scope: "evidence", src: source, mode: "list", title, rows: [], total: 0 };
      try { const p = new URLSearchParams({ source, cutoff: this.cutoff });
        const d = await getJSON("/api/evidence-aggregation/detail?" + p); if (req !== this._dreq) return;
        this.drill = { scope: "evidence", src: source, mode: "list", title, ...d };
      } catch (e) {}
    },
    closeDrill() { this.drill = null; this._dreq++; },

    // ----------------------------------------------- biological insights
    async loadInsights() {
      const req = ++this._nreq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { this.insights = null; return; }
        const d = await getJSON("/api/biological-insights?cutoff=" + this.cutoff);
        if (req !== this._nreq) return;
        this.insights = d;
      } catch (e) {}
    },
    openSection(sec) {
      this.section = sec; this.modalOpen = true; this.closeDrill();
      if (sec === "profiles") { if (!this.profiles) this.loadProfiles(); }
      else if (sec === "bio") { if (!this.bio) this.loadBio(); else this.$nextTick(() => this.drawSankey()); }
      else if (sec === "hubs") { if (!this.hubs) this.loadHubs(); }
      else if (sec === "impact") { if (!this.impact) this.loadImpact(); }
      else if (sec === "evidence") { if (!this.evid) this.loadEvidence(); }
      else if (sec === "enrich") { if (!this.enrich) this.loadEnrich(); }
      if (sec === "network") this.$nextTick(() => this.draw());
      if (sec === "bio") this.$nextTick(() => this.drawSankey());
    },
    closeSection() { this.modalOpen = false; },
    sectionTitle() {
      return ({ network: "Curated network", profiles: "Communication profiles", bio: "Biological processes",
                hubs: "Molecular hubs", impact: "Ligand–receptor rewiring", evidence: "Ligand–receptor dot plot",
                enrich: "Functional enrichment" })[this.section] || "";
    },
    async gotoInsight(ins) {
      const f = ins.focus || {};
      this.openSection(ins.goto);
      if (ins.goto === "profiles" && f.population) { if (!this.profiles) await this.loadProfiles(); this.profileDrill(f.population); }
      else if (ins.goto === "hubs" && f.name) { if (!this.hubs) await this.loadHubs(); this.hubDrill(f.kind, f.name); }
      else if (ins.goto === "bio" && f.source) { if (!this.bio) await this.loadBio(); this.bioSource = f.source; this.$nextTick(() => this.drawSankey()); if (f.name) this.selectProcess(f.name); }
    },

    // ----------------------------------------------- ligand–receptor dot plot
    async loadEvidence() {
      const req = ++this._ereq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/lr-dotplot?cutoff=" + this.cutoff);
        if (req !== this._ereq) return;
        this.evid = d;
        this.$nextTick(() => this.drawDotplot());
      } catch (e) {}
    },
    drawDotplot() {
      const d = this.evid; if (!d || !d.cells) return;
      const svg = document.getElementById("lr-dot"); if (!svg) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };
      const lrs = d.lr_order || [], cps = d.cp_order || [];
      const left = 250, top = 118, cw = 44, rh = 25;
      const W = left + cw * cps.length + 40, H = top + rh * lrs.length + 30;
      svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
      el("rect", { x: 0, y: 0, width: W, height: H, fill: "#ffffff" });
      const lrIx = {}; lrs.forEach((v, i) => lrIx[v] = i);
      const cpIx = {}; cps.forEach((v, i) => cpIx[v] = i);
      // column labels (rotated)
      cps.forEach((cp, j) => { const x = left + cw * j + cw / 2;
        el("text", { x, y: top - 8, transform: `rotate(-45 ${x} ${top - 8})`, "font-size": 10.5, fill: "#2a2d34" }).textContent = cp; });
      // row labels + gridlines
      lrs.forEach((lr, i) => { const y = top + rh * i + rh / 2;
        el("text", { x: left - 8, y: y + 4, "text-anchor": "end", "font-size": 11.5, "font-weight": 600, fill: "#2a2d34" }).textContent = lr;
        el("line", { x1: left, y1: y, x2: left + cw * cps.length, y2: y, stroke: "#eef0f3", "stroke-width": 1 });
      });
      // dots
      (d.cells || []).forEach(c => {
        const i = lrIx[c.lr], j = cpIx[c.cp]; if (i === undefined || j === undefined) return;
        const x = left + cw * j + cw / 2, y = top + rh * i + rh / 2;
        const rad = 3 + 5 * ((c.rank || 0) / 100);
        const rcp = (c.lr.split("→")[1] || c.lr);  // colour by receptor — method-independent
        el("circle", { cx: x, cy: y, r: rad.toFixed(1), fill: this._famColor(rcp), opacity: 0.85, stroke: "#fff", "stroke-width": 0.6 });
      });
    },
    // ----------------------------------------------- functional enrichment
    async loadEnrich() {
      const req = ++this._fereq; this.feLoading = true;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/functional-enrichment?cutoff=" + this.cutoff);
        if (req !== this._fereq) return;
        this.enrich = d;
        this.$nextTick(() => this.drawEnrich());
      } catch (e) { if (req === this._fereq) this.enrich = { sources: [], reached: false }; }
      finally { if (req === this._fereq) this.feLoading = false; }
    },
    _enrichBars(items, x0, y0, color, title, maxv, labw, barmax) {
      const F = 'font-family="Inter,Segoe UI,sans-serif"';
      let s = `<text x="${x0}" y="${y0 - 12}" ${F} font-size="13" font-weight="700" fill="#2a2d34">${title}</text>`;
      if (!items.length) { s += `<text x="${x0}" y="${y0 + 14}" ${F} font-size="11.5" fill="#9aa0ad">no significant terms</text>`; return s; }
      items.forEach((m, i) => {
        const y = y0 + i * 26;
        const nm = m.name.length > 40 ? m.name.slice(0, 39) + "…" : m.name;
        const bw = Math.max(3, barmax * m.neglog10p / (maxv || 1));
        s += `<text x="${x0 + labw - 6}" y="${y + 13}" text-anchor="end" ${F} font-size="11" fill="#2a2d34">${nm.replace(/&/g, "&amp;").replace(/</g, "&lt;")}</text>`;
        s += `<rect x="${x0 + labw}" y="${y + 3}" width="${bw.toFixed(1)}" height="15" rx="3" fill="${color}"/>`;
        s += `<text x="${(x0 + labw + bw + 5).toFixed(1)}" y="${y + 14}" ${F} font-size="10.5" fill="#6b6f78">${m.neglog10p}</text>`;
      });
      return s;
    },
    drawEnrich() {
      const d = this.enrich; if (!d || !d.sources) return;
      d.sources.forEach((src, idx) => {
        const svg = document.getElementById("fe-" + idx); if (!svg) return;
        while (svg.firstChild) svg.removeChild(svg.firstChild);
        const rows = Math.max(src.before.length, src.after.length, 1);
        const W = 1180, colW = 560, labw = 300, barmax = 190, y0 = 34;
        const H = y0 + rows * 26 + 20;
        svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
        const maxv = Math.max(1, ...src.before.map(x => x.neglog10p), ...src.after.map(x => x.neglog10p));
        const inner = this._enrichBars(src.before, 10, y0, "#9aa0ad", "BEFORE · all candidate genes", maxv, labw, barmax)
          + this._enrichBars(src.after, 10 + colW + 40, y0, "#e4572e", "AFTER · curated genes", maxv, labw, barmax);
        svg.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>${inner}`;
      });
    },
    evPct() { return this.evid && this.evid.total_kept ? Math.round((this.evid.computed / this.evid.total_kept) * 100) : 0; },
    evBar(v, total) { return Math.round((v / (total || 1)) * 100) + "%"; },
    srcLabel(k) { const s = (this.evid ? this.evid.sources : []).find(x => x.key === k); return s ? s.label : k; },
    cellCls(v) { return v === true ? "yes" : v === false ? "no" : "na"; },
    cellMark(v) { return v === true ? "✓" : v === false ? "·" : "–"; },
    drawUpset() {
      const svg = document.getElementById("ev-upset");
      if (!svg || !this.evid) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const combos = (this.evid.combos || []).slice(0, 14);
      const keys = (this.evid.sources || []).map(s => s.key);
      if (!combos.length || !keys.length) return;
      const el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };
      const L = 150, colW = 48, rowH = 26, barH = 120, pad = 14, dotR = 5;
      const W = L + combos.length * colW + pad, H = barH + 24 + keys.length * rowH + 10;
      svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
      const maxC = Math.max(...combos.map(c => c.count), 1);
      const cx = i => L + colW * i + colW / 2;
      const matrixTop = barH + 24;

      // bars
      combos.forEach((c, i) => {
        const h = Math.max(2, (c.count / maxC) * (barH - 18));
        el("rect", { x: cx(i) - 12, y: barH - h, width: 24, height: h, rx: 2, fill: "var(--layer-data)" });
        const t = el("text", { x: cx(i), y: barH - h - 5, "text-anchor": "middle", class: "ev-barnum" });
        t.textContent = c.count;
      });
      // source row labels + guide stripes
      keys.forEach((k, r) => {
        const y = matrixTop + rowH * r + rowH / 2;
        if (r % 2 === 0) el("rect", { x: L - 4, y: y - rowH / 2, width: W - L, height: rowH, fill: "var(--bg-inset)", opacity: "0.5" });
        const t = el("text", { x: L - 14, y: y + 4, "text-anchor": "end", class: "ev-rowlbl" });
        t.textContent = this.srcLabel(k);
      });
      // dot matrix
      combos.forEach((c, i) => {
        const on = new Set(c.sources);
        const ys = keys.map((k, r) => ({ k, y: matrixTop + rowH * r + rowH / 2, on: on.has(k) }));
        const lit = ys.filter(d => d.on);
        if (lit.length > 1) el("line", { x1: cx(i), y1: lit[0].y, x2: cx(i), y2: lit[lit.length - 1].y, stroke: "var(--text-strong)", "stroke-width": 2 });
        ys.forEach(d => el("circle", { cx: cx(i), cy: d.y, r: dotR, fill: d.on ? "var(--text-strong)" : "var(--border)" }));
        if (!lit.length) { const t = el("text", { x: cx(i), y: matrixTop - 4, "text-anchor": "middle", class: "ev-none" }); t.textContent = "none"; }
      });
    },
    pctKept(metric) {
      if (!this.impact || !this.impact.before) return "";
      const b = this.impact.before[metric] || 0, a = this.impact.after[metric] || 0;
      return b ? Math.round((a / b) * 100) + "%" : "–";
    },
    drawMiniNet(svgId, net) {
      const svg = document.getElementById(svgId);
      if (!svg || !net) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const W = 420, H = 420, cx = W / 2, cy = H / 2, R = 138;
      const nodes = net.nodes.filter(d => d.out + d.in > 0);
      const n = nodes.length || 1;
      const pos = {};
      nodes.forEach((d, i) => { const a = -Math.PI / 2 + (2 * Math.PI * i) / n; pos[d.name] = { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a), a }; });
      const el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };
      const defs = el("defs", {});
      const m = el("marker", { id: svgId + "-arr", viewBox: "0 0 10 10", refX: "9", refY: "5", markerWidth: "6", markerHeight: "6", markerUnits: "userSpaceOnUse", orient: "auto-start-reverse" }, defs);
      el("path", { d: "M 0 0 L 10 5 L 0 10 z", fill: "var(--layer-data)" }, m);
      const edges = net.edges.slice(0, 120);
      const maxC = Math.max(...edges.map(e => e.count), 1);
      const maxDeg = Math.max(...nodes.map(d => d.out + d.in), 1);
      const rOf = {}; nodes.forEach(d => rOf[d.name] = 4 + 9 * Math.sqrt((d.out + d.in) / maxDeg));
      const g = el("g", {});
      edges.forEach(e => {
        const s = pos[e.sender], t = pos[e.receiver]; if (!s || !t || e.sender === e.receiver) return;
        const dx = t.x - s.x, dy = t.y - s.y, len = Math.hypot(dx, dy) || 1, ux = dx / len, uy = dy / len;
        const sx = s.x + ux * rOf[e.sender], sy = s.y + uy * rOf[e.sender], tx = t.x - ux * (rOf[e.receiver] + 2), ty = t.y - uy * (rOf[e.receiver] + 2);
        const mx = (sx + tx) / 2, my = (sy + ty) / 2, bend = 0.18 * len;
        el("path", { d: `M ${sx} ${sy} Q ${mx - uy * bend} ${my + ux * bend} ${tx} ${ty}`, fill: "none", stroke: "var(--layer-data)",
          "stroke-width": 0.5 + 3 * (e.count / maxC), "stroke-opacity": 0.22 + 0.5 * (e.count / maxC), "marker-end": `url(#${svgId}-arr)` }, g);
      });
      nodes.forEach(d => { const s = pos[d.name];
        el("circle", { cx: s.x, cy: s.y, r: rOf[d.name], fill: "var(--bg-shell)", stroke: "var(--text-strong)", "stroke-width": 1 }, g);
        const ca = Math.cos(s.a), sa = Math.sin(s.a);
        const lx = cx + (R + rOf[d.name] + 9) * ca, ly = cy + (R + rOf[d.name] + 9) * sa;
        const anchor = Math.abs(ca) < 0.3 ? "middle" : (ca > 0 ? "start" : "end");
        const tx = el("text", { x: lx, y: ly + 4, "text-anchor": anchor, class: "mini-label" }, g);
        tx.textContent = this._trunc(d.name, 16);
        const tt = el("title", {}, g); tt.textContent = `${d.name}: sends ${d.out}, receives ${d.in}`;
      });
    },

    async load() {
      const req = ++this._req;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { this.data = null; return; }
        const d = await getJSON("/api/curated-network?" + this.params());
        if (req !== this._req) return;
        this.data = d;
        this._counts = {};
        d.edges.forEach(e => { this._counts[e.sender + "|" + e.receiver] = e.count; });
        this.$nextTick(() => this.draw());
        if (this.sel) this.openPair(this.sel.sender, this.sel.receiver);
      } catch (e) {}
    },

    shownTypes() { return this.data ? this.data.nodes.map(n => n.name) : []; },
    isShown(t) { return !this.f.types.length || this.f.types.includes(t); },
    toggleType(t) {
      let cur = this.f.types.length ? [...this.f.types] : [...this.data.all_cell_types];
      cur = cur.includes(t) ? cur.filter(x => x !== t) : [...cur, t];
      this.f.types = cur.length === this.data.all_cell_types.length ? [] : cur;
      this.load();
    },
    count(s, r) { return this._counts[s + "|" + r] || 0; },
    maxCount() { return this.data && this.data.edges.length ? this.data.edges[0].count : 1; },
    cellStyle(s, r) {
      const c = this.count(s, r);
      if (!c) return "";
      const a = 0.12 + 0.88 * Math.sqrt(c / this.maxCount());
      return `background: color-mix(in srgb, var(--layer-data) ${Math.round(a * 100)}%, transparent); color: ${a > 0.55 ? "#fff" : "var(--text-strong)"}`;
    },

    async openPair(s, r) {
      this.sel = { sender: s, receiver: r }; this.pairData = null;
      const p = new URLSearchParams({ sender: s, receiver: r, cutoff: this.cutoff });
      try { this.pairData = await getJSON("/api/curated-network/pair?" + p); } catch (e) {}
      this.draw();
    },

    draw() {
      const svg = document.getElementById("nw-svg");
      if (!svg || !this.data) return;
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const W = 760, H = 720, cx = W / 2, cy = H / 2, R = 250;
      const nodes = this.data.nodes, n = nodes.length;
      const pos = {};
      nodes.forEach((d, i) => {
        const a = -Math.PI / 2 + (2 * Math.PI * i) / Math.max(n, 1);
        pos[d.name] = { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a), a };
      });
      const el = (tag, attrs, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent.appendChild(e); return e; };

      const defs = el("defs", {});
      for (const [id, col] of [["arr", "var(--layer-data)"], ["arr-out", "var(--accent)"], ["arr-in", "var(--layer-data)"],
                               ["arr-sel", "var(--text-strong)"], ["arr-dim", "var(--border)"]]) {
        const m = el("marker", { id, viewBox: "0 0 10 10", refX: "9", refY: "5", markerWidth: "7", markerHeight: "7",
                                  markerUnits: "userSpaceOnUse", orient: "auto-start-reverse" }, defs);
        el("path", { d: "M 0 0 L 10 5 L 0 10 z", fill: col }, m);
      }

      const edges = this.data.edges.slice(0, this.maxEdges);
      const maxC = Math.max(...edges.map(e => e.count), 1);
      const maxDeg = Math.max(...nodes.map(d => d.out + d.in), 1);
      const nodeR = d => 7 + 13 * Math.sqrt((d.out + d.in) / maxDeg);
      const rOf = {}; nodes.forEach(d => rOf[d.name] = nodeR(d));

      const gEdges = el("g", {});
      const edgeEls = [];
      edges.forEach(e => {
        const s = pos[e.sender], t = pos[e.receiver];
        if (!s || !t) return;
        const w = 0.6 + 5.5 * (e.count / maxC);
        let d;
        if (e.sender === e.receiver) {
          const ox = Math.cos(s.a), oy = Math.sin(s.a), r0 = rOf[e.sender];
          const p1x = s.x + ox * r0, p1y = s.y + oy * r0;
          d = `M ${p1x - oy * 6} ${p1y + ox * 6} C ${p1x + ox * 38 - oy * 26} ${p1y + oy * 38 + ox * 26}, ${p1x + ox * 38 + oy * 26} ${p1y + oy * 38 - ox * 26}, ${p1x + oy * 6} ${p1y - ox * 6}`;
        } else {
          const dx = t.x - s.x, dy = t.y - s.y, len = Math.hypot(dx, dy) || 1;
          const ux = dx / len, uy = dy / len;
          const sx = s.x + ux * rOf[e.sender], sy = s.y + uy * rOf[e.sender];
          const tx = t.x - ux * (rOf[e.receiver] + 3), ty = t.y - uy * (rOf[e.receiver] + 3);
          const mx = (sx + tx) / 2, my = (sy + ty) / 2;
          const bend = 0.18 * len;
          d = `M ${sx} ${sy} Q ${mx - uy * bend} ${my + ux * bend} ${tx} ${ty}`;
        }
        const isSel = this.sel && this.sel.sender === e.sender && this.sel.receiver === e.receiver;
        const p = el("path", { d, fill: "none", stroke: isSel ? "var(--text-strong)" : "var(--layer-data)",
          "stroke-width": w, "stroke-opacity": isSel ? 1 : 0.28 + 0.5 * (e.count / maxC),
          "marker-end": `url(#${isSel ? "arr-sel" : "arr"})`, class: "nw-edge" }, gEdges);
        const tip = el("title", {}, p); tip.textContent = `${e.sender} → ${e.receiver}: ${e.count} kept interactions`;
        p.addEventListener("click", () => this.openPair(e.sender, e.receiver));
        edgeEls.push({ p, e, w, isSel });
      });

      const gNodes = el("g", {});
      nodes.forEach(d => {
        const s = pos[d.name];
        const g = el("g", { class: "nw-node", tabindex: "0" }, gNodes);
        el("circle", { cx: s.x, cy: s.y, r: rOf[d.name], fill: "var(--bg-shell)", stroke: "var(--text-strong)", "stroke-width": 1.5 }, g);
        const lx = cx + (R + rOf[d.name] + 10) * Math.cos(s.a), ly = cy + (R + rOf[d.name] + 10) * Math.sin(s.a) + 8 * Math.sin(s.a);
        const right = Math.cos(s.a) >= -0.01;
        const txt = el("text", { x: lx, y: ly + 4, "text-anchor": Math.abs(Math.cos(s.a)) < 0.12 ? "middle" : (right ? "start" : "end"), class: "nw-label" }, g);
        txt.textContent = d.name;
        const t2 = el("title", {}, g); t2.textContent = `${d.name}: sends ${d.out}, receives ${d.in} (${d.self} within the population)`;
        const hi = on => edgeEls.forEach(({ p, e, w, isSel }) => {
          const out = e.sender === d.name, inc = e.receiver === d.name;
          if (!on) { p.setAttribute("stroke", isSel ? "var(--text-strong)" : "var(--layer-data)"); p.setAttribute("stroke-opacity", isSel ? 1 : 0.28 + 0.5 * (e.count / maxC)); p.setAttribute("marker-end", `url(#${isSel ? "arr-sel" : "arr"})`); return; }
          p.setAttribute("stroke", out ? "var(--accent)" : inc ? "var(--layer-data)" : "var(--border)");
          p.setAttribute("stroke-opacity", out || inc ? 0.9 : 0.15);
          p.setAttribute("marker-end", `url(#${out ? "arr-out" : inc ? "arr-in" : "arr-dim"})`);
        });
        g.addEventListener("mouseenter", () => hi(true)); g.addEventListener("mouseleave", () => hi(false));
        g.addEventListener("focus", () => hi(true)); g.addEventListener("blur", () => hi(false));
      });
    },

    // ------------------------------------------------ figure export (SVG + PNG)
    // Exports each view for assembly in draw.io. The four native views
    // (network, before/after networks, sankey) are serialised straight from
    // their live <svg> with all computed styles inlined, so they open as clean
    // editable vectors. The HTML views (hubs, processes, profiles) are redrawn
    // as vector SVG from the same data the page shows.
    _palette() {
      const cs = getComputedStyle(document.documentElement);
      const g = (v, d) => { const x = (cs.getPropertyValue(v) || "").trim(); return x || d; };
      return { accent: g("--accent", "#e4572e"), data: g("--layer-data", "#3f7cbf"),
               text: g("--text-strong", "#1f2328"), muted: "#6b6f78" };
    },
    _dl(blob, name) {
      const u = URL.createObjectURL(blob);
      const a = document.createElement("a"); a.href = u; a.download = name;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(u), 2000);
    },
    _sizeOf(svg) {
      const vb = (svg.getAttribute("viewBox") || "0 0 800 600").split(/\s+/).map(Number);
      return { w: vb[2] || 800, h: vb[3] || 600 };
    },
    _pngFromSvg(svgStr, w, h, name, scale = 3) {
      const img = new Image();
      img.onload = () => {
        const c = document.createElement("canvas"); c.width = Math.round(w * scale); c.height = Math.round(h * scale);
        const ctx = c.getContext("2d"); ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, c.width, c.height);
        ctx.setTransform(scale, 0, 0, scale, 0, 0); ctx.drawImage(img, 0, 0, w, h);
        c.toBlob(b => { if (b) this._dl(b, name + ".png"); }, "image/png");
      };
      img.onerror = () => {};
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svgStr);
    },
    _serializeInlineSvg(svgId) {
      const src = document.getElementById(svgId);
      if (!src || !src.childNodes.length) return null;
      const clone = src.cloneNode(true);
      const props = ["fill", "fill-opacity", "stroke", "stroke-width", "stroke-opacity",
                     "stroke-dasharray", "stroke-linecap", "stroke-linejoin", "opacity",
                     "font-family", "font-size", "font-weight", "font-style",
                     "text-anchor", "dominant-baseline", "marker-start", "marker-end", "marker-mid"];
      const walk = (a, b) => {
        const cs = getComputedStyle(a);
        let s = "";
        for (const p of props) { const v = cs.getPropertyValue(p); if (v) s += p + ":" + v + ";"; }
        b.setAttribute("style", s);
        const ac = a.children, bc = b.children;
        for (let i = 0; i < ac.length; i++) if (bc[i]) walk(ac[i], bc[i]);
      };
      walk(src, clone);
      const { w, h } = this._sizeOf(src);
      clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
      clone.setAttribute("width", w); clone.setAttribute("height", h);
      const bg = document.createElementNS(NS, "rect");
      bg.setAttribute("x", 0); bg.setAttribute("y", 0); bg.setAttribute("width", w); bg.setAttribute("height", h);
      bg.setAttribute("fill", "#ffffff");
      clone.insertBefore(bg, clone.firstChild);
      const raw = new XMLSerializer().serializeToString(clone);
      return { str: '<?xml version="1.0" encoding="UTF-8"?>\n' + this._resolveVars(raw), w, h };
    },
    // Resolve any leftover CSS custom properties (e.g. inside <defs>/markers that
    // getComputedStyle skips because they are not rendered) to concrete values.
    _resolveVars(str) {
      const cs = getComputedStyle(document.documentElement);
      return str.replace(/var\(\s*(--[\w-]+)\s*(?:,[^)]*)?\)/g, (m, name) => {
        const v = (cs.getPropertyValue(name) || "").trim();
        return v || m;
      });
    },
    exportSVG(svgId, name) {
      const r = this._serializeInlineSvg(svgId);
      if (!r) { alert("Nothing to export yet — open this view first."); return; }
      this._dl(new Blob([r.str], { type: "image/svg+xml" }), name + ".svg");
      this._pngFromSvg(r.str, r.w, r.h, name);
    },
    _esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); },
    _wrapSvg(inner, w, h) {
      return '<?xml version="1.0" encoding="UTF-8"?>\n' +
        `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}">` +
        `<rect x="0" y="0" width="${w}" height="${h}" fill="#ffffff"/>${inner}</svg>`;
    },
    _niceNum(x) {
      if (!(x > 0)) return 1;
      const e = Math.floor(Math.log10(x)), f = x / Math.pow(10, e);
      const nf = f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10;
      return nf * Math.pow(10, e);
    },
    _niceTicks(max, n) {
      const step = this._niceNum(max / n);
      const out = [];
      for (let v = 0; v <= max + 1e-9; v += step) out.push(Math.round(v * 10) / 10);
      if (out.length < 2) out.push(Math.round(max * 10) / 10);
      return out;
    },
    _barColumn(items, maxv, x0, y0, color, title, labw, barmax, P, rowh = 26, axis = "") {
      const F = 'font-family="Inter,Segoe UI,sans-serif"';
      const plotX = x0 + labw, rows = Math.max(items.length, 1);
      const plotTop = y0 - 6, plotBot = y0 + rows * rowh - 6;
      let s = `<text x="${x0}" y="${y0 - 16}" ${F} font-size="11" font-weight="700" letter-spacing="0.06em" fill="${P.muted}">${this._esc(title).toUpperCase()}</text>`;
      this._niceTicks(maxv, 4).forEach(t => {
        const gx = plotX + barmax * t / (maxv || 1);
        s += `<line x1="${gx.toFixed(1)}" y1="${plotTop}" x2="${gx.toFixed(1)}" y2="${plotBot}" stroke="#eceef1" stroke-width="1"/>`;
        s += `<text x="${gx.toFixed(1)}" y="${plotBot + 15}" text-anchor="middle" ${F} font-size="9.5" fill="#aab0b8">${t}</text>`;
      });
      if (axis) s += `<text x="${plotX + barmax / 2}" y="${plotBot + 31}" text-anchor="middle" ${F} font-size="10" fill="${P.muted}">${this._esc(axis)}</text>`;
      items.forEach((m, i) => {
        const y = y0 + i * rowh;
        const val = m.interactions != null ? m.interactions : m.ccis;
        const bw = Math.max(2, barmax * val / (maxv || 1));
        const nm = m.name.length > 30 ? m.name.slice(0, 29) + "…" : m.name;
        s += `<text x="${plotX - 8}" y="${y + 13}" text-anchor="end" ${F} font-size="11.5" fill="#3a3f47">${this._esc(nm)}</text>`;
        s += `<rect x="${plotX}" y="${y + 3}" width="${bw.toFixed(1)}" height="14" rx="3" fill="${color}"/>`;
        s += `<text x="${(plotX + bw + 6).toFixed(1)}" y="${y + 13}" ${F} font-size="10.5" fill="${P.muted}">${val}</text>`;
      });
      return s;
    },
    exportHubs(name) {
      if (!this.hubs || (!this.hubs.ligands.length && !this.hubs.receptors.length)) { alert("No hub data to export yet."); return; }
      const P = this._palette();
      const L = this.hubs.ligands.slice(0, 12), R = this.hubs.receptors.slice(0, 12);
      const rows = Math.max(L.length, R.length), w = 1060, y0 = 54, h = y0 + rows * 26 + 40;
      let inner = this._barColumn(L, this.hubs.max_ligand, 30, y0, "#3E7C8A", "Top ligands", 150, 300, P, 26, "kept interactions");
      inner += this._barColumn(R, this.hubs.max_receptor, 560, y0, "#5A7D9A", "Top receptors", 150, 300, P, 26, "kept interactions");
      const svg = this._wrapSvg(inner, w, h);
      this._dl(new Blob([svg], { type: "image/svg+xml" }), name + ".svg");
      this._pngFromSvg(svg, w, h, name);
    },
    exportProc(name) {
      const b = this.bioCur ? this.bioCur() : null;
      if (!b || !b.processes || !b.processes.length) { alert("No process data to export yet."); return; }
      const P = this._palette();
      const items = b.processes.slice(0, 14), w = 680, y0 = 54, h = y0 + items.length * 26 + 40;
      const inner = this._barColumn(items, b.max, 30, y0, "#6A8E8A", "Processes by kept interactions", 220, 360, P, 26, "kept interactions");
      const svg = this._wrapSvg(inner, w, h);
      this._dl(new Blob([svg], { type: "image/svg+xml" }), name + ".svg");
      this._pngFromSvg(svg, w, h, name);
    },
    exportProfiles(name) {
      if (!this.profiles || !this.profiles.profiles.length) { alert("No profile data to export yet."); return; }
      const P = this._palette(), F = 'font-family="Inter,Segoe UI,sans-serif"';
      const ps = this.profiles.profiles, rowh = 40, w = 680, h = 80 + ps.length * rowh;
      let inner = `<text x="30" y="34" ${F} font-size="15" font-weight="700" fill="${P.text}">Communication profiles — send / receive</text>`;
      const barx = 180, barmax = 360, maxTotal = Math.max.apply(null, ps.map(p => p.total || 1));
      ps.forEach((p, i) => {
        const y = 56 + i * rowh, tot = p.total || (p.outgoing + p.incoming) || 1;
        const bw = Math.max(6, barmax * tot / maxTotal);
        const outW = bw * (p.outgoing / ((p.outgoing + p.incoming) || 1));
        inner += `<text x="${barx - 10}" y="${y + 15}" text-anchor="end" ${F} font-size="12.5" font-weight="600" fill="${P.text}">${this._esc(p.name)}</text>`;
        inner += `<rect x="${barx}" y="${y + 3}" width="${outW.toFixed(1)}" height="16" fill="${P.accent}"/>`;
        inner += `<rect x="${(barx + outW).toFixed(1)}" y="${y + 3}" width="${(bw - outW).toFixed(1)}" height="16" fill="${P.data}"/>`;
        inner += `<text x="${(barx + bw + 6).toFixed(1)}" y="${y + 15}" ${F} font-size="11" fill="${P.muted}">${p.outgoing} send · ${p.incoming} recv</text>`;
      });
      inner += `<rect x="30" y="${h - 24}" width="12" height="12" fill="${P.accent}"/><text x="48" y="${h - 14}" ${F} font-size="11" fill="${P.muted}">sends</text>`;
      inner += `<rect x="112" y="${h - 24}" width="12" height="12" fill="${P.data}"/><text x="130" y="${h - 14}" ${F} font-size="11" fill="${P.muted}">receives</text>`;
      const svg = this._wrapSvg(inner, w, h);
      this._dl(new Blob([svg], { type: "image/svg+xml" }), name + ".svg");
      this._pngFromSvg(svg, w, h, name);
    },
  });
}
