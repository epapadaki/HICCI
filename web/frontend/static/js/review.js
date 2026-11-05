/* Review queue page. The queue is computed on the server; opening an event
   shows the same documentation panel as the Interactions page, starting on
   the Assessment tab. */

function reviewApp() {
  const DEFAULTS = { single: 99.9 };

  return Object.assign(docPanel(), {
    data: { categories: [], methods: [], rows: [], total: 0, page: 1, pages: 1 },
    category: null,
    page: 1,
    t: { ...DEFAULTS },
    criteriaOpen: false,
    sortBy: "",
    loading: false,
    selectedEvent: null,
    _req: 0,

    async init() {
      const p = new URLSearchParams(location.search);
      this.category = p.get("category");
      this.sortBy = p.get("sort_by") || "";
      this.page = parseInt(p.get("page") || "1", 10);
      for (const k of Object.keys(DEFAULTS)) if (p.get(k) !== null) this.t[k] = parseFloat(p.get(k));
      this.initPanel();
      this.loadAssessed();
      await this.load();
      const sel = p.get("selected");
      if (sel !== null) {
        this.selectedEvent = p.get("event");
        this.select(parseInt(sel, 10), { keepUrl: true });
      }
    },

    params() {
      const p = new URLSearchParams();
      if (this.category) p.set("category", this.category);
      if (this.sortBy) p.set("sort_by", this.sortBy);
      if (this.page > 1) p.set("page", this.page);
      for (const [k, v] of Object.entries(this.t)) if (v !== DEFAULTS[k]) p.set(k, v);
      return p;
    },

    writeUrl() {
      const p = this.params();
      if (this.selectedId !== null) { p.set("selected", this.selectedId); if (this.selectedEvent) p.set("event", this.selectedEvent); }
      const qs = p.toString();
      history.replaceState(null, "", qs ? "?" + qs : location.pathname);
    },

    async load() {
      const req = ++this._req;
      this.loading = true;
      const p = this.params();
      p.set("page", this.page);
      try {
        const d = await getJSON("/api/review-queue?" + p.toString());
        if (req !== this._req) return;
        this.data = d;
        this.page = d.page;
      } catch (e) {
        if (req === this._req) this.data = { ...this.data, rows: [], total: 0 };
      } finally {
        if (req === this._req) this.loading = false;
      }
      this.writeUrl();
    },

    onAssessmentSaved() { this.load(); },

    open(row, tab) {
      this.selectedEvent = row.event_key;
      this.openDoc(row.result_id, tab || "assessment");
    },

    setCategory(c) { this.category = c; this.page = 1; this.load(); },
    goPage(n) { this.page = n; this.load(); const w = document.querySelector(".table-wrap"); if (w) w.scrollTop = 0; },
    resetCriteria() { this.t = { ...DEFAULTS }; this.page = 1; this.load(); },
    totalQueued() { return (this.data.categories || []).reduce((s, c) => s + c.count, 0); },
    catLabel(key) { const c = (this.data.categories || []).find(x => x.key === key); return c ? c.label : key; },
    closePanel() { this.selectedId = null; this.selectedEvent = null; this.writeUrl(); },
  });
}
