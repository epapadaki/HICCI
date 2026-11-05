/* Shared by all pages: the biological documentation panel, the expert
   assessment form, the experiment description and small helpers.
   Each page merges docPanel() into its own Alpine component and provides
   writeUrl() (and optionally onAssessmentSaved()). */

async function getJSON(url, options) {
  const res = await fetch(url, options);
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch (e) {}
    const err = new Error(detail);
    err.status = res.status;
    throw err;
  }
  return res.json();
}

function docPanel() {
  return {
    theme: document.documentElement.getAttribute("data-theme") || "light",
    selectedId: null,
    tab: "biology", bioSub: "mechanism",
    doc: { overview: null, biology: null, biologyError: null, literature: null, literatureError: null, assessments: null, evidence: null, evidenceError: null, mech: null, response: null, responseError: null },
    _docReq: 0,
    showAll: {},

    assessed: {},
    keepCutoff: 95,   // default: pre-select interactions with rank-in-tool >= this
    reasons: {},
    form: { expert: "", decision: "", relevance: 3, confidence: 3, reasons: [], comment: "" },
    formMsg: { text: "", error: false },
    saving: false,

    datasetsInfo: { active: null, items: [] },
    resultsState: null, tools: [],
    docOpen: false,
    DOC_TABS: [
      { key: "biology", icon: "ti-dna", label: "Prior knowledge", short: "Prior knowledge" },
      { key: "literature", icon: "ti-book", label: "Literature support", short: "Literature support" },
      { key: "relevance", icon: "ti-target-arrow", label: "Case study relevance", short: "Case study relevance" },
      { key: "data", icon: "ti-chart-bar", label: "Case study", short: "Case study" },
      { key: "response", icon: "ti-activity-heartbeat", label: "Signalling response", short: "Signalling" },
    ],
    mechSignal: true, mechBig: false, mechView: "schematic", mechRx: 0, rxFailed: false,

    context: {},
    contextDraft: {},
    contextOpen: false,

    initPanel() {
      this.loadContext();
      getJSON("/api/datasets").then(d => { this.datasetsInfo = d; }).catch(() => {});
      this.watchResults();
      try { this.form.expert = localStorage.getItem("ccc-expert") || ""; } catch (e) {}
      try { const c = parseFloat(localStorage.getItem("ccc-keep-cutoff")); if (!isNaN(c)) this.keepCutoff = c; } catch (e) {}
    },
    async watchResults() {
      // Poll until the results index is ready (or known to be missing), for the top-bar notice.
      let wasBuilding = false;
      for (let i = 0; i < 400; i++) {
        try { const st = await getJSON("/api/status"); this.resultsState = st.state; if (st.methods) this.tools = st.methods; } catch (e) { return; }
        if (this.resultsState !== "building") {
          // The page asked for data while the index was being built: load it again now that it is ready.
          if (wasBuilding && this.resultsState === "ready" && !location.pathname.startsWith("/datasets")) location.reload();
          return;
        }
        wasBuilding = true;
        await new Promise(r => setTimeout(r, 1500));
      }
    },

    writeUrl() {},

    async switchDataset(id) {
      if (!id || id === this.datasetsInfo.active) return;
      await getJSON("/api/datasets/active", {
        method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ id }),
      });
      // Start clean (filters and selection belong to the old dataset); from the
      // Datasets page, continue to the interactions of the chosen dataset.
      location.href = location.pathname === "/datasets" ? "/interactions" : location.pathname;
    },
    activeDataset() { return (this.datasetsInfo.items || []).find(d => d.id === this.datasetsInfo.active) || null; },
    async unloadDataset() {
      await getJSON("/api/datasets/clear", { method: "PUT" });
      if (location.pathname === "/datasets") location.reload(); else location.href = "/datasets";
    },
    onAssessmentSaved() {},

    // ------------------------------------------------------------ documentation
    select(id, opts = {}) {
      const req = ++this._docReq;
      this.selectedId = id;
      this.doc = { overview: null, biology: null, biologyError: null, literature: null, literatureError: null, assessments: null, evidence: null, evidenceError: null, mech: null, response: null, responseError: null };
      this.resetForm();
      this.showAll = {};
      this.bioSub = "mechanism";
      if (!opts.keepUrl) this.writeUrl();
      const panel = document.getElementById("doc-panel");
      if (panel) panel.scrollTop = 0;

      getJSON(`/api/results/${id}`)
        .then(d => { if (req === this._docReq) this.doc.overview = d; })
        .catch(e => { if (req === this._docReq) this.doc.biologyError = "Could not load this result: " + e.message; });

      getJSON(`/api/results/${id}/biology`)
        .then(d => { if (req === this._docReq) this.doc.biology = d; })
        .catch(e => { if (req === this._docReq) this.doc.biologyError = "Could not load the biological annotations: " + e.message; });

      this.loadLiterature(false, req);
      if (this.tab === "response") this.$nextTick(() => this.loadResponse());
      this.loadAssessments(req);

      getJSON(`/api/results/${id}/mechanism`)
        .then(d => { if (req === this._docReq) { this.doc.mech = d; this.mechRx = 0; this.rxFailed = false; } })
        .catch(() => {});

      getJSON(`/api/results/${id}/dataset-evidence`)
        .then(d => { if (req === this._docReq) this.doc.evidence = d; })
        .catch(e => { if (req === this._docReq) this.doc.evidenceError = "Could not load expression data: " + e.message; });
    },

    // ------------------------------------------------------------ expert assessment (on/off)
    loadAssessed() {
      getJSON("/api/assessments/keys").then(d => { this.assessed = d; }).catch(() => {});
    },
    // Best within-tool rank for a row (ranks are not comparable across tools, so
    // we use the row's strongest tool rank). 0–100; higher = stronger.
    rowBestRank(row) {
      if (Array.isArray(row.tools)) { const rs = row.tools.map(t => t.rank).filter(r => r != null); return rs.length ? Math.max(...rs) : 0; }
      return (row.global_percentile == null ? 0 : row.global_percentile);
    },
    // Default: interactions ranked >= keepCutoff in their tool are kept; an
    // explicit choice (the expert ticked or unticked the box) always wins.
    isKept(row) {
      const d = this.assessed[row.event_key];
      if (d === "Accept") return true;
      if (d === "Reject") return false;
      return this.rowBestRank(row) >= this.keepCutoff;
    },
    rememberReviewer() { try { localStorage.setItem("ccc-expert", (this.form.expert || "").trim()); } catch (e) {} },
    setKeepCutoff(v) {
      let c = parseFloat(v); if (isNaN(c)) c = 95; c = Math.max(0, Math.min(100, c));
      this.keepCutoff = c;
      try { localStorage.setItem("ccc-keep-cutoff", String(c)); } catch (e) {}
    },
    async toggleKept(row, checked) {
      const id = (row.open_id != null) ? row.open_id : row.id;
      this.assessed[row.event_key] = checked ? "Accept" : "Reject";   // store an explicit override
      this.assessed = { ...this.assessed };
      try {
        await getJSON(`/api/results/${id}/keep`, {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ kept: checked, expert: this.form.expert }),
        });
        try { if ((this.form.expert || "").trim()) localStorage.setItem("ccc-expert", this.form.expert.trim()); } catch (e) {}
      } catch (e) {
        alert("Could not save: " + e.message);
        this.loadAssessed();   // revert to the server's truth
      }
    },
    // Delete stored decisions so the rank cutoff governs pre-selection again.
    // kind: null = all, "Accept" or "Reject" = only that kind.
    async clearDecisions(kind) {
      const what = kind === "Accept" ? "accepted" : kind === "Reject" ? "rejected" : "all";
      if (!confirm(`Clear ${what} expert decisions for this dataset?\n\nThe rank cutoff will govern selection again. This cannot be undone.`)) return;
      try {
        const r = await getJSON("/api/assessments/clear", {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ decision: kind || null }),
        });
        this.assessed = r.keys || {};
      } catch (e) { alert("Could not clear: " + e.message); }
    },
    loadAssessments(req = this._docReq) {
      getJSON(`/api/results/${this.selectedId}/assessments`)
        .then(d => { if (req === this._docReq) { this.doc.assessments = d; this.reasons = d.reasons; } })
        .catch(() => {});
    },
    resetForm() {
      this.form = { expert: this.form.expert, decision: "", relevance: 3, confidence: 3, reasons: [], comment: "" };
      this.formMsg = { text: "", error: false };
    },
    async saveAssessment() {
      if (!this.form.decision) { this.formMsg = { text: "Choose a decision first.", error: true }; return; }
      if (!this.form.expert.trim()) { this.formMsg = { text: "Enter your name or identifier.", error: true }; return; }
      this.saving = true;
      try {
        const res = await getJSON(`/api/results/${this.selectedId}/assessments`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(this.form),
        });
        try { localStorage.setItem("ccc-expert", this.form.expert.trim()); } catch (e) {}
        this.resetForm();
        this.formMsg = { text: `Saved (assessment ${res.id}).`, error: false };
        this.loadAssessments();
        this.loadAssessed();
        this.onAssessmentSaved();
      } catch (e) {
        this.formMsg = { text: "Not saved: " + e.message, error: true };
      } finally {
        this.saving = false;
      }
    },
    pct(v) { return v === null || v === undefined ? "–" : (v * 100 >= 99.95 ? "100" : (v * 100).toFixed(v * 100 < 10 ? 1 : 0)) + "%"; },
    ordinal(n) { const s = ["th", "st", "nd", "rd"], v = n % 100; return n + (s[(v - 20) % 10] || s[v] || s[0]); },

    loadResponse() {
      if (this.doc.response || this.selectedId === null || this._respLoading === this.selectedId) return;
      const id = this.selectedId, req = this._docReq;
      this._respLoading = id;
      getJSON(`/api/results/${id}/response`)
        .then(d => { if (req === this._docReq) this.doc.response = d; })
        .catch(e => { if (req === this._docReq) this.doc.responseError = e.message; })
        .finally(() => { if (this._respLoading === id) this._respLoading = null; });
    },
    barStyle(e) { const w = Math.min(50, Math.abs(e) / 3 * 50); return (e >= 0 ? "left:50%;" : "right:50%;") + "width:" + w + "%"; },
    fmtEffect(p) { return (p.effect_size >= 0 ? "+" : "") + p.effect_size.toFixed(2) + (p.fdr < 0.05 ? " *" : ""); },
    fmtFdr(q) { return q < 0.001 ? "<0.001" : q.toFixed(3); },

    activeRuns() {
      const d = this.activeDataset();
      return d && d.runs ? d.runs.filter(r => r.available) : [];
    },
    async switchRun(run) {
      try {
        await getJSON(`/api/datasets/${encodeURIComponent(this.datasetsInfo.active)}/view`, {
          method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ run }) });
        location.href = location.pathname;
      } catch (e) { alert(e.message); }
    },

    shortTool(m) { return { CellPhoneDB: "CPDB", CellChat: "CellChat", NATMI: "NATMI", SingleCellSignalR: "SCSR", ICELLNET: "ICELLNET" }[m] || m; },
    nTools() { return this.tools.length || 3; },
    toolsRange() { return Array.from({ length: this.nTools() }, (_, i) => i + 1); },

    fmtDate(s) { return (s || "").replace("T", " ").slice(0, 16); },

    loadLiterature(refresh = false, req = this._docReq) {
      const id = this.selectedId;
      if (id === null) return;
      if (refresh) { this.doc.literature = null; this.doc.literatureError = null; }
      getJSON(`/api/results/${id}/literature${refresh ? "?refresh=true" : ""}`)
        .then(d => { if (req === this._docReq) this.doc.literature = d; })
        .catch(e => { if (req === this._docReq) this.doc.literatureError = "Could not load the literature: " + e.message; });
    },

    openDoc(id, tab) {
      if (id !== this.selectedId || !this.doc.overview) this.select(id);
      this.tab = tab || this.tab || "biology";
      this.docOpen = true;
      if (this.tab === "response") this.$nextTick(() => this.loadResponse());
      this.$nextTick(() => document.querySelectorAll(".doc-grid > [role=tabpanel]").forEach(b => { b.scrollTop = 0; }));
    },
    setTab(t) { this.tab = t; if (t === "response") this.loadResponse(); document.querySelectorAll(".doc-grid > [role=tabpanel]").forEach(b => { b.scrollTop = 0; }); },
    closeDoc() {
      this.docOpen = false;
      const u = new URL(location.href);
      if (u.searchParams.has("selected")) { u.searchParams.delete("selected"); history.replaceState(null, "", u); }
    },
    closePanel() { this.closeDoc(); },
    onEscape() { if (this.contextOpen) this.contextOpen = false; else if (this.mechBig) this.mechBig = false; else this.closeDoc(); },

    entityNote(entity, side) {
      if (!entity) return "";
      const who = side === "ligand" ? "Ligand" : "Receptor";
      if (entity.kind === "complex") return `${who} is a complex of ${entity.subunits.join(", ")}; all subunits are kept.`;
      if (entity.kind === "proxy") return `${who} is a small molecule; the method estimates it from the expression of the enzyme that makes it.`;
      if (entity.kind === "other") return `${who} is a CellPhoneDB entity without verified subunits.`;
      return "";
    },

    // ------------------------------------------------------------ experiment context
    async loadContext() {
      try { this.context = await getJSON("/api/project-context"); } catch (e) { this.context = {}; }
    },
    contextFilled() {
      const c = this.context || {};
      return Boolean(c.tissue_or_organ || c.disease_or_condition || (c.experimental_groups || []).length);
    },
    contextLabel() {
      const c = this.context || {};
      const parts = [c.organism, c.tissue_or_organ].filter(Boolean);
      return parts.length ? parts.join(", ") : "Not described";
    },
    contextSummary() {
      const c = this.context || {};
      return [c.organism, c.tissue_or_organ, c.disease_or_condition, ...(c.experimental_groups || [])].filter(Boolean).join(", ");
    },
    openContext() {
      const c = this.context || {};
      this.contextDraft = {
        organism: c.organism || "",
        tissue_or_organ: c.tissue_or_organ || "",
        disease_or_condition: c.disease_or_condition || "",
        groupsText: (c.experimental_groups || []).join(", "),
        notes: c.notes || "",
      };
      this.contextOpen = true;
      this.$nextTick(() => { const el = document.getElementById("c-org"); if (el) el.focus(); });
    },
    async saveContext() {
      const d = this.contextDraft;
      this.context = await getJSON("/api/project-context", {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          organism: d.organism,
          tissue_or_organ: d.tissue_or_organ,
          disease_or_condition: d.disease_or_condition,
          experimental_groups: d.groupsText,
          notes: d.notes,
        }),
      });
      this.contextOpen = false;
      if (this.selectedId !== null) { this.doc.literature = null; this.loadLiterature(); }
    },

    // ------------------------------------------------------------ misc
    toggleTheme() {
      this.theme = this.theme === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", this.theme);
      try { localStorage.setItem("ccc-theme", this.theme); } catch (e) {}
    },

    groupPathways(list) {
      const groups = new Map();
      (list || []).forEach(pw => {
        const key = pw.name.toLowerCase();
        if (!groups.has(key)) groups.set(key, { name: pw.name, items: [], directions: [] });
        const g = groups.get(key);
        if (!g.items.some(i => i.pmid === pw.pmid)) g.items.push(pw);
        if (!g.directions.includes(pw.direction)) g.directions.push(pw.direction);
      });
      return [...groups.values()].sort((a, b) => b.items.length - a.items.length);
    },

    dataVerdict() {
      const e = this.doc.evidence;
      if (!e || !e.available) return null;
      const sideDet = (s) => {
        if (s.joint) return s.joint.detected;
        const f = s.subunits.filter(x => x.found).map(x => x.detected);
        return f.length ? Math.min(...f) : 0;
      };
      const lig = sideDet(e.ligand), rcp = sideDet(e.receptor);
      const found = e.receptor.subunits.filter(x => x.found);
      let rank = null, nPop = null;
      if (found.length) { const lim = found.reduce((a, b) => (a.rank >= b.rank ? a : b)); rank = lim.rank; nPop = lim.n_populations; }
      const topFrac = (rank && nPop) ? rank / nPop : null;
      const present = lig >= 0.1 && rcp >= 0.1;
      const specific = topFrac !== null && topFrac <= 0.34;
      let level = "ok", text;
      if (!present) { level = "weak"; text = "weakly present in your data"; }
      else if (specific) { level = "strong"; text = "present and fairly specific to these populations"; }
      else { level = "ok"; text = "present, but broadly expressed across populations"; }
      return { lig, rcp, rank, nPop, present, specific, level, text };
    },

    relevanceSummary() {
      const pubs = (this.doc.literature && this.doc.literature.publications) || [];
      const n = pubs.length;
      const pairSpecific = pubs.filter(p => p.pair_specific).length;
      const coMention = pubs.filter(p => p.co_mention_only).length;
      // how many studies match the user's context on each aspect
      const aspectMatch = {};
      pubs.forEach(p => (p.relevance && p.relevance.items || []).forEach(a => {
        if (a.status === "same" || a.status === "partly") {
          aspectMatch[a.aspect] = (aspectMatch[a.aspect] || 0) + 1;
        }
      }));
      // a study "close to your experiment": pair-specific AND at least one context match
      const close = pubs.filter(p => p.pair_specific &&
        (p.relevance && p.relevance.items || []).some(a => a.status === "same" || a.status === "partly")).length;
      return { n, pairSpecific, coMention, aspectMatch, close };
    },

    relevanceRank(p) {
      // pair-specific first, then those matching the user's context, newest first
      const grp = p.pair_specific ? 0 : (p.co_mention_only ? 2 : 1);
      const ctx = (p.relevance && p.relevance.items || []).some(a => a.status === "same" || a.status === "partly") ? 0 : 1;
      const yr = parseInt(p.year, 10); return [grp, ctx, isNaN(yr) ? 0 : -yr];
    },

    relevanceSorted() {
      const pubs = ((this.doc.literature && this.doc.literature.publications) || []).slice();
      return pubs.sort((a, b) => { const x = this.relevanceRank(a), y = this.relevanceRank(b);
        for (let i = 0; i < 3; i++) if (x[i] !== y[i]) return x[i] - y[i]; return 0; });
    },

    statusLabel(s) {
      return {
        same: "Same", partly: "Partly", different: "Different", no_name_match: "No name match",
        not_reported: "Not in the study", not_set: "Not in your description", mentioned: "Mentioned",
        not_mentioned: "Not mentioned", info: "Reported",
      }[s] || s;
    },

    fmtInt(n) { return (n ?? 0).toLocaleString("en-US"); },
    fmtPct(v) { return v === null || v === undefined ? "–" : (v >= 99.95 ? "100" : Number(v).toFixed(1)); },
    fmtScore(v) {
      if (v === null || v === undefined) return "–";
      const a = Math.abs(v);
      if (a !== 0 && (a < 0.001 || a >= 10000)) return Number(v).toExponential(2);
      return Number(v).toPrecision(3).replace(/\.?0+$/, "");
    },
    // DOMParser does not run scripts or load images, so this is safe for
    // text coming from local databases.
    stripHtml(s) { return new DOMParser().parseFromString(s || "", "text/html").body.textContent || ""; },
  };
}
