/* Interactions page: the results table. The table asks the server for one
   page at a time; the documentation panel (common.js) loads its sections
   independently. */

function cccApp() {
  const DEFAULT_FILTERS = { q: "", sender: "", receiver: "", methods: [], minMethods: 1, minPercentile: 0 };

  return Object.assign(docPanel(), {
    index: { state: "idle", message: "" },
    facets: { methods: [], senders: [], receivers: [] },
    filters: { ...DEFAULT_FILTERS, methods: [] },
    sort: "tools",
    view: "interactions",          // "interactions": one row per interaction; "tools": one row per tool prediction
    selectedEvent: null,
    results: { total: 0, page: 1, pages: 1, rows: [], score_sort_allowed: false },
    loading: false,
    loadingMore: false,
    _resultsReq: 0,

    // ------------------------------------------------------------ startup
    async init() {
      this.readUrl();
      this.initPanel();
      await this.waitForIndex();
      if (this.index.state !== "ready") return;   // dataset without results: the banner explains
      this.loadAssessed();
      await this.loadFacets();
      await this.loadResults();
      if (this.selectedId !== null) { this.select(this.selectedId, { keepUrl: true }); this.docOpen = true; }
      window.addEventListener("popstate", () => { this.readUrl(); this.applyFilters(); });
      this.setupInfiniteScroll();
    },

    async waitForIndex() {
      for (;;) {
        try { this.index = await getJSON("/api/status"); } catch (e) { this.index = { state: "error", message: e.message }; }
        if (this.index.state !== "building") return;   // ready, missing, no_dataset or error
        await new Promise(r => setTimeout(r, 1000));
      }
    },

    async loadFacets() {
      try { this.facets = await getJSON("/api/facets"); } catch (e) { /* banner shows state */ }
    },

    // ------------------------------------------------------------ url state
    readUrl() {
      const p = new URLSearchParams(location.search);
      this.filters = {
        q: p.get("q") || "",
        sender: p.get("sender") || "",
        receiver: p.get("receiver") || "",
        methods: p.getAll("method"),
        minMethods: parseInt(p.get("min_methods") || "1", 10),
        minPercentile: parseFloat(p.get("min_percentile") || "0"),
      };
      this.view = p.get("view") === "tools" ? "tools" : "interactions";
      this.sort = p.get("sort") || (this.view === "tools" ? "percentile" : "tools");
      this.results.page = 1;
      const sel = p.get("selected");
      this.selectedId = sel !== null ? parseInt(sel, 10) : null;
    },

    queryParams() {
      const p = new URLSearchParams();
      const f = this.filters;
      if (f.q) p.set("q", f.q);
      if (f.sender) p.set("sender", f.sender);
      if (f.receiver) p.set("receiver", f.receiver);
      f.methods.forEach(m => p.append("method", m));
      if (f.minMethods > 1) p.set("min_methods", f.minMethods);
      if (f.minPercentile > 0) p.set("min_percentile", f.minPercentile);
      if (this.view === "tools") p.set("view", "tools");
      if (this.sort !== (this.view === "tools" ? "percentile" : "tools")) p.set("sort", this.sort);
      return p;
    },

    writeUrl() {
      const p = this.queryParams();
      if (this.selectedId !== null) p.set("selected", this.selectedId);
      const qs = p.toString();
      history.replaceState(null, "", qs ? "?" + qs : location.pathname);
    },

    // ------------------------------------------------------------ results
    async loadResults(append = false) {
      const req = ++this._resultsReq;
      if (append) this.loadingMore = true; else this.loading = true;
      try {
        const p = this.queryParams();
        p.set("page", this.results.page);
        const url = this.view === "tools" ? "/api/results?" : "/api/interactions?";
        const data = await getJSON(url + p.toString());
        if (req !== this._resultsReq) return;           // a newer request superseded this one
        if (append) { this.results = { ...data, rows: this.results.rows.concat(data.rows) }; }
        else { this.results = data; }
        this.sort = data.sort;
      } catch (e) {
        if (req === this._resultsReq && !append) this.results = { total: 0, page: 1, pages: 1, rows: [], score_sort_allowed: false };
      } finally {
        if (req === this._resultsReq) { this.loading = false; this.loadingMore = false; }
      }
      this.writeUrl();
    },

    // Load the next page and append it (infinite scroll — no page buttons).
    loadMore() {
      if (this.loading || this.loadingMore) return;
      if (!this.results.pages || this.results.page >= this.results.pages) return;
      this.results.page += 1;
      this.loadResults(true);
    },
    setupInfiniteScroll() {
      if (this._onScroll) return;
      const wrap = document.querySelector(".table-wrap");
      const near = () => {
        if (wrap && wrap.scrollHeight > wrap.clientHeight + 4) {
          return wrap.scrollTop + wrap.clientHeight >= wrap.scrollHeight - 800;   // table-wrap scrolls
        }
        const doc = document.documentElement;                                      // otherwise the page scrolls
        return window.innerHeight + window.scrollY >= doc.scrollHeight - 800;
      };
      this._onScroll = () => { if (near()) this.loadMore(); };
      if (wrap) wrap.addEventListener("scroll", this._onScroll, { passive: true });
      window.addEventListener("scroll", this._onScroll, { passive: true });
    },

    applyFilters() {
      this.results.page = 1;
      const w = document.querySelector(".table-wrap"); if (w) w.scrollTop = 0; window.scrollTo(0, 0);
      this.loadResults();
    },
    toggleMethod(m) {
      const i = this.filters.methods.indexOf(m);
      if (i >= 0) this.filters.methods.splice(i, 1); else this.filters.methods.push(m);
      this.applyFilters();
    },
    hasFilters() {
      const f = this.filters;
      return f.q || f.sender || f.receiver || f.methods.length || f.minMethods > 1 || f.minPercentile > 0;
    },
    resetFilters() { this.filters = { ...DEFAULT_FILTERS, methods: [] }; this.sort = "percentile"; this.applyFilters(); },
    // Pathway-level curation: Accept / Reject every interaction matching the current filters.
    bulkBusy: false,
    async bulkAssess(kept) {
      const total = this.results.total || 0;
      if (!total) { alert("No interactions match the current filters."); return; }
      const verb = kept ? "Accept" : "Reject";
      const f = this.filters;
      const body = { expert: this.form.expert, filters: {
        methods: f.methods, sender: f.sender, receiver: f.receiver,
        q: f.q, min_percentile: f.minPercentile, min_methods: f.minMethods } };
      this.bulkBusy = true;
      try {
        // 1) preview: show exactly what will change, broken down by signalling category
        const pv = await getJSON("/api/assessments/bulk", {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ ...body, kept, dry_run: true }),
        });
        const lines = (pv.by_category || []).map(([c, n]) => `   • ${c}: ${n}`).join("\n");
        if (!confirm(`${verb} ${pv.total} interaction(s)?\n\nBy signalling category:\n${lines}\n\nThese are exactly the rows currently listed in the table.`)) {
          this.bulkBusy = false; return;
        }
        // 2) apply
        const r = await getJSON("/api/assessments/bulk", {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ ...body, kept, dry_run: false }),
        });
        try { if ((this.form.expert || "").trim()) localStorage.setItem("ccc-expert", this.form.expert.trim()); } catch (e) {}
        await this.loadAssessed();
        alert(`${verb}ed ${r.count} interaction(s).` + (r.errors ? ` (${r.errors} could not be saved.)` : ""));
      } catch (e) {
        alert("Bulk action failed: " + e.message);
      } finally {
        this.bulkBusy = false;
      }
    },
    setSort(key) { this.sort = key; this.applyFilters(); },
    // ---- pathway-level curation panel ----
    pathways: [],
    pathwaysOpen: false,
    pathBusy: false,
    pathQ: "",
    async togglePathwayPanel() {
      this.pathwaysOpen = !this.pathwaysOpen;
      if (this.pathwaysOpen) await this.loadPathways();
    },
    async loadPathways() {
      this.pathBusy = true;
      try { const d = await getJSON("/api/pathways"); this.pathways = d.pathways || []; }
      catch (e) { alert("Could not load pathways: " + e.message); }
      finally { this.pathBusy = false; }
    },
    pathStatus(p) {
      if (p.kept && !p.rejected) return "kept";
      if (p.rejected && !p.kept) return "rejected";
      if (p.kept && p.rejected) return "mixed";
      return "undecided";
    },
    pathList() {
      const q = (this.pathQ || "").toLowerCase().trim();
      return q ? this.pathways.filter(p => p.category.toLowerCase().includes(q)) : this.pathways;
    },
    keptPathwayCount() { return this.pathways.filter(p => this.pathStatus(p) === "kept").length; },
    keptInteractionCount() { return this.pathways.reduce((s, p) => s + (p.kept || 0), 0); },
    async decidePathway(p, kept) {
      this.pathBusy = true;
      try {
        await getJSON("/api/pathways/decide", {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ category: p.category, kept, expert: this.form.expert }),
        });
        try { if ((this.form.expert || "").trim()) localStorage.setItem("ccc-expert", this.form.expert.trim()); } catch (e) {}
        await this.loadPathways();
        await this.loadAssessed();
        this.loadResults();
      } catch (e) { alert("Failed: " + e.message); }
      finally { this.pathBusy = false; }
    },
    setView(v) {
      if (this.view === v) return;
      this.view = v; this.sort = v === "tools" ? "percentile" : "tools";
      this.results = { total: 0, page: 1, pages: 1, rows: [], score_sort_allowed: false };
      this.applyFilters();
    },
    openEvent(row, id, tab) { this.selectedEvent = row.event_key; this.openDoc(id || row.open_id, tab || "biology"); },
    toolRank(row, m) { const t = row.tools.find(x => x.method === m); return t ? t : { rank: null, id: null }; },
    sortState(key) { return this.sort === key ? (["sender", "receiver", "ligand", "receptor"].includes(key) ? "ascending" : "descending") : "none"; },
    goPage(n) {
      this.results.page = Math.max(1, Math.min(n, this.results.pages));
      this.loadResults();
      const wrap = document.querySelector(".table-wrap");
      if (wrap) wrap.scrollTop = 0;
    },

  });
}
