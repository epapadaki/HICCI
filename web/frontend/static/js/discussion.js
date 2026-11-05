/* Discussion page: biological insights + LLM-written, PubMed-grounded synthesis
   of the curated results and their functional enrichment. */
function cccDiscussion() {
  return Object.assign(docPanel(), {
    insights: null,
    _nreq: 0,
    cutoff: 95,
    dStatus: null,
    disc: null,
    dLoading: false,
    dCopied: false,

    async init() {
      this.initPanel();
      try { this.cutoff = parseInt(localStorage.getItem("ccc-keep-cutoff") || "95", 10) || 95; } catch (e) {}
      this.loadInsights();
      this.loadDStatus();
    },

    setCutoff(v) {
      v = Math.max(0, Math.min(100, parseInt(v, 10) || 0));
      this.cutoff = v;
      try { localStorage.setItem("ccc-keep-cutoff", String(v)); } catch (e) {}
      this.loadInsights();
    },

    async loadInsights() {
      const req = ++this._nreq;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { if (req === this._nreq) this.insights = null; return; }
        const d = await getJSON("/api/biological-insights?cutoff=" + this.cutoff);
        if (req === this._nreq) this.insights = d;
      } catch (e) { if (req === this._nreq) this.insights = null; }
    },

    async loadDStatus() {
      try {
        this.dStatus = await getJSON("/api/discussion/status");
      } catch (e) { this.dStatus = { available: false }; }
    },

    async genDiscussion(refresh) {
      this.dLoading = true; this.dCopied = false;
      try {
        const st = await getJSON("/api/status");
        if (st.state !== "ready") { this.disc = { note: "No prepared results yet." }; return; }
        this.disc = await getJSON("/api/discussion?cutoff=" + this.cutoff + (refresh ? "&refresh=1" : ""));
      } catch (e) {
        this.disc = { error: "Request failed: " + (e.message || e) };
      } finally { this.dLoading = false; }
    },

    copyDisc() {
      if (!this.disc || !this.disc.text) return;
      navigator.clipboard.writeText(this.disc.text).then(() => {
        this.dCopied = true; setTimeout(() => this.dCopied = false, 1500);
      }).catch(() => {});
    },

    renderDisc(src) {
      const esc = s => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
      const cite = s => s.replace(/\[(\d+(?:\s*,\s*\d+)*)\]/g, (m, grp) =>
        "<sup>[" + grp.split(",").map(x => x.trim()).map(n => `<a href="#ref-${n}">${n}</a>`).join(", ") + "]</sup>");
      const inline = s => cite(esc(s).replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>"));
      return String(src || "").trim().split(/\n{2,}/).map(b => {
        const t = b.trim();
        if (/^#{1,3}\s+/.test(t)) return `<h3>${inline(t.replace(/^#{1,3}\s+/, ""))}</h3>`;
        return `<p>${inline(t).replace(/\n/g, "<br>")}</p>`;
      }).join("");
    },
  });
}
