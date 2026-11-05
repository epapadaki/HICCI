/* Method comparison page. */
function methodsApp() {
  return Object.assign(docPanel(), {
    s: null, error: null, combo: null, ev: null,

    async init() {
      this.initPanel();
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { this.error = st.message || "The results of this dataset are not ready yet."; return; }
        this.s = await getJSON("/api/method-comparison");
        const all = this.s.combinations.find(c => c.methods.length === this.s.methods.length && c.count > 0)
                 || this.s.combinations.find(c => c.count > 0);
        if (all) { this.combo = all.key; this.loadEvents(1); }     // the list is never empty
      } catch (e) { this.error = e.message; }
    },

    short(m) { return { CellPhoneDB: "CPDB", CellChat: "CChat", SingleCellSignalR: "SCSR", ICELLNET: "ICELL" }[m] || m; },
    pctOf(v, total) { return total ? Math.max(0.5, v / total * 100) : 0; },
    maxCombo() { return Math.max(...this.s.combinations.map(c => c.count), 1); },

    async pickCombo(key) {
      this.combo = key; this.ev = null;
      await this.loadEvents(1);
      const el = document.getElementById("combo-events"); if (el) el.scrollIntoView({ behavior: "smooth", block: "start" });
    },

    async loadEvents(page) {
      try { this.ev = await getJSON(`/api/method-comparison/events?combo=${encodeURIComponent(this.combo)}&page=${page}`); } catch (e) {}
    },
  });
}
