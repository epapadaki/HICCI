/* Receptor explorer page. */
function receptorsApp() {
  const PAGE = 50;
  return Object.assign(docPanel(), {
    data: null, f: { top: "all", min_events: 5 }, q: "", sort: "fdr", page: 1,
    methodsList: ["CellPhoneDB", "NATMI", "CellChat"],
    detail: null, tab: "ligands", fn: null, fnBusy: false, fnError: null, _req: 0,

    async init() { this.initPanel(); this.$watch("q", () => { this.page = 1; }); await this.load(); },
    params() { return new URLSearchParams({ top: this.f.top, min_events: this.f.min_events }); },
    async load() {
      const req = ++this._req;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") return;
        const d = await getJSON("/api/receptors?" + this.params());
        if (req !== this._req) return;
        if (d.methods && d.methods.length === 1 && (this.f.top === "all" || this.f.top === "two")) {
          this.methodsList = d.methods; this.f.top = d.methods[0] + ":95"; return this.load();   // one tool: its own top
        }
        this.data = d; this.page = 1; if (d.methods) this.methodsList = d.methods;
        if (this.detail) this.open(this.detail.receptor_key);
        else if (d.receptors.length) this.open(d.receptors[0].receptor_key);
      } catch (e) {}
    },
    filtered() {
      if (!this.data) return [];
      const q = this.q.trim().toUpperCase();
      let rows = q ? this.data.receptors.filter(r => r.receptor.toUpperCase().includes(q)) : [...this.data.receptors];
      const key = this.sort;
      rows.sort((a, b) => key === "fdr" ? (a.fdr - b.fdr) || ((b.fold || 0) - (a.fold || 0)) : ((b[key] || 0) - (a[key] || 0)));
      return rows;
    },
    pages() { return Math.max(1, Math.ceil(this.filtered().length / PAGE)); },
    pageRows() { return this.filtered().slice((this.page - 1) * PAGE, this.page * PAGE); },
    sortBy(k) { this.sort = k; this.page = 1; },
    foldWidth(f) { if (!f) return 0; return Math.min(100, Math.log2(Math.max(f, 1)) / Math.log2(32) * 100); },
    fmtP(p) { if (p === null || p === undefined) return "–"; if (p < 1e-4) return "<0.0001"; return p < 0.01 ? p.toFixed(4) : p.toFixed(3); },
    async open(key) {
      this.fn = null; this.fnError = null;
      if (!this.detail || this.detail.receptor_key !== key) this.tab = "ligands";
      try { this.detail = await getJSON(`/api/receptors/detail?key=${encodeURIComponent(key)}&` + this.params()); } catch (e) {}
    },
    async loadFunctions() {
      this.fnBusy = true; this.fnError = null;
      try { this.fn = await getJSON(`/api/receptors/ligand-enrichment?key=${encodeURIComponent(this.detail.receptor_key)}&` + this.params()); }
      catch (e) { this.fnError = e.message; }
      this.fnBusy = false;
    },
  });
}
