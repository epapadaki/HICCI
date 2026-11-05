/* Functional enrichment page: before vs after over-representation, local. */
function cccEnrich() {
  const NS = "http://www.w3.org/2000/svg";
  return Object.assign(docPanel(), {
    enrich: null,
    feLoading: false,
    _fereq: 0,
    cutoff: 95,
    showShift: false,

    async init() {
      this.initPanel();
      try { this.cutoff = parseInt(localStorage.getItem("ccc-keep-cutoff") || "95", 10) || 95; } catch (e) {}
      this.loadEnrich();
      window.addEventListener("resize", () => this.drawEnrich());
    },

    setCutoff(v) {
      v = Math.max(0, Math.min(100, parseInt(v, 10) || 0));
      this.cutoff = v;
      try { localStorage.setItem("ccc-keep-cutoff", String(v)); } catch (e) {}
      this.loadEnrich();
    },

    async loadEnrich() {
      const req = ++this._fereq;
      this.feLoading = true;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { if (req === this._fereq) this.enrich = { sources: [], reached: false, note: "No prepared results for this dataset yet." }; return; }
        const d = await getJSON("/api/functional-enrichment?cutoff=" + this.cutoff + (this.showShift ? "&descriptive=1" : ""));
        if (req !== this._fereq) return;
        this.enrich = d;
        this.$nextTick(() => this.drawEnrich());
      } catch (e) {
        if (req === this._fereq) this.enrich = { sources: [], reached: false, note: "Could not compute enrichment." };
      } finally {
        if (req === this._fereq) this.feLoading = false;
      }
    },

    toggleShift() {
      this.showShift = !this.showShift;
      if (this.showShift && this.enrich && !this.enrich.descriptive) this.loadEnrich();
      else this.$nextTick(() => this.drawEnrich());
    },

    emptyNote() {
      if (!this.enrich || !this.enrich.sources) return "";
      const noTerms = this.enrich.sources.filter(s => s.installed !== false && (!s.after || !s.after.length)).map(s => s.label);
      const notInst = this.enrich.sources.filter(s => s.installed === false).map(s => s.label);
      const parts = [];
      if (noTerms.length) parts.push("No specifically-enriched terms in: " + noTerms.join(", ") + ".");
      if (notInst.length) parts.push("Not installed: " + notInst.join(", ") + ".");
      return parts.join("  ");
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

    _enrichBars(items, x0, y0, color, title, maxv, labw, barmax, rows) {
      const F = 'font-family="Inter,Segoe UI,sans-serif"';
      const esc = t => String(t).replace(/&/g, "&amp;").replace(/</g, "&lt;");
      const plotX = x0 + labw, plotTop = y0 - 6, plotBot = y0 + rows * 26 - 6;
      let s = `<text x="${x0}" y="${y0 - 16}" ${F} font-size="11" font-weight="700" letter-spacing="0.06em" fill="#8a8f98">${esc(title.toUpperCase())}</text>`;
      // gridlines + tick labels + axis title (behind the bars)
      this._niceTicks(maxv, 4).forEach(t => {
        const gx = plotX + barmax * t / (maxv || 1);
        s += `<line x1="${gx.toFixed(1)}" y1="${plotTop}" x2="${gx.toFixed(1)}" y2="${plotBot}" stroke="#eceef1" stroke-width="1"/>`;
        s += `<text x="${gx.toFixed(1)}" y="${plotBot + 15}" text-anchor="middle" ${F} font-size="9.5" fill="#aab0b8">${t}</text>`;
      });
      s += `<text x="${plotX + barmax / 2}" y="${plotBot + 31}" text-anchor="middle" ${F} font-size="10" fill="#8a8f98">−log10(FDR)</text>`;
      if (!items.length) { s += `<text x="${x0}" y="${y0 + 14}" ${F} font-size="11.5" fill="#aab0b8">no significant terms</text>`; return s; }
      items.forEach((m, i) => {
        const y = y0 + i * 26;
        const nm = m.name.length > 44 ? m.name.slice(0, 43) + "…" : m.name;
        const bw = Math.max(2, barmax * m.neglog10p / (maxv || 1));
        s += `<text x="${plotX - 8}" y="${y + 13}" text-anchor="end" ${F} font-size="11" fill="#3a3f47">${esc(nm)}</text>`;
        s += `<rect x="${plotX}" y="${y + 3}" width="${bw.toFixed(1)}" height="14" rx="3" fill="${color}"/>`;
        s += `<text x="${(plotX + bw + 6).toFixed(1)}" y="${y + 13.5}" ${F} font-size="10" fill="#8a8f98">${m.neglog10p}</text>`;
      });
      return s;
    },

    drawEnrich() {
      const d = this.enrich;
      if (!d || !d.sources) return;
      const W = 1180, y0 = 38;
      const AFTER = "#3E7C8A", BEFORE = "#b6bcc6";
      d.sources.forEach((src, idx) => {
        // Inferential panel (curated vs candidate background)
        const svg = document.getElementById("fe-" + idx);
        if (svg) {
          if (src.installed === false) {
            const H = 56;
            svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
            svg.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>`
              + `<text x="10" y="32" font-family="Inter,Segoe UI,sans-serif" font-size="12" fill="#9aa0ad">Library not installed — run  python scripts/fetch_genesets.py  then reload.</text>`;
          } else {
            const rows = Math.max(src.after.length, 1);
            const H = y0 + rows * 26 + 44;
            svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
            const labw = 470, barmax = 600;
            const maxv = Math.max(1, ...src.after.map(x => x.neglog10p));
            const inner = this._enrichBars(src.after, 10, y0, AFTER, "Curated genes · enriched vs all candidate L–R genes", maxv, labw, barmax, rows);
            svg.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>${inner}`;
          }
        }
        // Descriptive 'compositional shift' panel (before vs after, genome background)
        const svg2 = document.getElementById("fs-" + idx);
        if (svg2 && src.before_genome !== undefined) {
          const bg = src.before_genome || [], ag = src.after_genome || [];
          const rows = Math.max(bg.length, ag.length, 1);
          const H = y0 + rows * 26 + 44;
          svg2.setAttribute("viewBox", `0 0 ${W} ${H}`);
          const colW = 560, labw = 300, barmax = 190;
          const maxv = Math.max(1, ...bg.map(x => x.neglog10p), ...ag.map(x => x.neglog10p));
          const inner = this._enrichBars(bg, 10, y0, BEFORE, "Default · all candidate genes", maxv, labw, barmax, rows)
            + this._enrichBars(ag, 10 + colW + 40, y0, AFTER, "Curated · expert-reviewed genes", maxv, labw, barmax, rows);
          svg2.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>${inner}`;
        }
      });
    },

    _dl(blob, name) {
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = name;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    },

    exportSVG(svgId, name) {
      const svg = document.getElementById(svgId);
      if (!svg || !svg.childNodes.length) { alert("Nothing to export yet."); return; }
      const vb = (svg.getAttribute("viewBox") || "0 0 1180 400").split(/\s+/).map(Number);
      const w = vb[2] || 1180, h = vb[3] || 400;
      const clone = svg.cloneNode(true);
      // Strip Alpine/framework binding attributes (:id, :aria-label, x-…, @…) and
      // inline style — a leftover ":" attribute is invalid XML and silently breaks
      // the SVG→PNG rasterisation (only the .svg would download otherwise).
      Array.prototype.slice.call(clone.attributes).forEach(a => {
        if (a.name[0] === ":" || a.name[0] === "@" || a.name.indexOf("x-") === 0 || a.name === "style") clone.removeAttribute(a.name);
      });
      clone.setAttribute("xmlns", NS);
      clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
      clone.setAttribute("viewBox", "0 0 " + w + " " + h);
      clone.setAttribute("width", w);
      clone.setAttribute("height", h);
      const raw = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
      this._dl(new Blob([raw], { type: "image/svg+xml" }), name + ".svg");
      const img = new Image();
      img.onload = () => {
        const c = document.createElement("canvas");
        const sc = 3; c.width = Math.round(w * sc); c.height = Math.round(h * sc);
        const ctx = c.getContext("2d");
        ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, c.width, c.height);
        ctx.setTransform(sc, 0, 0, sc, 0, 0); ctx.drawImage(img, 0, 0, w, h);
        c.toBlob(b => { if (b) this._dl(b, name + ".png"); }, "image/png");
      };
      img.onerror = () => {};
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(raw);
    },
  });
}
