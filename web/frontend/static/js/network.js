/* Network page: circular layout of cell populations, drawn as SVG. */
function networkApp() {
  const NS = "http://www.w3.org/2000/svg";
  return Object.assign(docPanel(), {
    data: null, view: "network", maxEdges: 80, typesOpen: false,
    f: { source: "consensus", min_methods: 3, min_rank: 90, types: [] },
    sel: null, pairData: null, _req: 0, _counts: {},

    async init() { this.initPanel(); await this.load(); window.addEventListener("resize", () => this.draw()); },

    params() {
      const p = new URLSearchParams({ source: this.f.source, min_methods: this.f.min_methods, min_rank: this.f.min_rank });
      this.f.types.forEach(t => p.append("cell_type", t));
      return p;
    },

    async load() {
      const req = ++this._req;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { this.data = null; return; }
        const d = await getJSON("/api/network?" + this.params());
        if (req !== this._req) return;
        if (d.methods && d.methods.length === 1 && this.f.source === "consensus") {
          this.f.source = d.methods[0]; return this.load();          // one tool: show its own results
        }
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
      const p = new URLSearchParams({ sender: s, receiver: r, source: this.f.source, min_methods: this.f.min_methods, min_rank: this.f.min_rank });
      try { this.pairData = await getJSON("/api/network/pair?" + p); } catch (e) {}
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
          const bend = 0.18 * len;                       // curve to one side: A→B and B→A separate
          d = `M ${sx} ${sy} Q ${mx - uy * bend} ${my + ux * bend} ${tx} ${ty}`;
        }
        const isSel = this.sel && this.sel.sender === e.sender && this.sel.receiver === e.receiver;
        const p = el("path", { d, fill: "none", stroke: isSel ? "var(--text-strong)" : "var(--layer-data)",
          "stroke-width": w, "stroke-opacity": isSel ? 1 : 0.28 + 0.5 * (e.count / maxC),
          "marker-end": `url(#${isSel ? "arr-sel" : "arr"})`, class: "nw-edge" }, gEdges);
        const tip = el("title", {}, p); tip.textContent = `${e.sender} → ${e.receiver}: ${e.count} interactions`;
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
  });
}
