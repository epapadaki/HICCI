/* Datasets page: analysis environment, dataset list, add-a-dataset steps. */

function datasetsApp() {
  const EMPTY = { h5ad: "", metadata: "", matrix: "X", label: "", id: "", organism: "", tissue: "",
                  cell_type_column: "", sample_column: "", group_column: "" };

  return Object.assign(docPanel(), {
    list: [],
    details: {},
    env: { checked: false, ok: false, message: "", path: "" },
    envEdit: false, envDraft: "", envSaving: false,
    step: 1, busy: false, msg: { text: "", error: false },
    form: { ...EMPTY }, idTouched: false, info: {},
    job: {}, lastId: null,
    methods: null, rDraft: "", cpdbDraft: "", cpdbSaving: false, runJob: null, checking: {}, checkingAll: false,
    ncbiEmail: "", ncbiKey: "", ncbiKeySet: false, ncbiSaved: "",
    ephemeral: true, ephemeralSaved: "",
    flow: null, setupOpen: false, files: { h5ad: [], annotation: [] },
    h5adName: "", h5adSize: null, metaName: "", dragOver: false, dragOverMeta: false, uploading: false, upPct: 0,
    toolList: ["CellPhoneDB", "NATMI", "CellChat", "SingleCellSignalR", "ICELLNET"], runTools: [],
    advice: {}, runMode: "single", singleTool: "", rec: null,
    TOOL_INFO: [
      { name: "CellPhoneDB", icon: "ti-phone-call", logo: "/static/img/tools/cellphonedb.png", tagline: "Curated database with multi-subunit complexes; fast.", features: ["Complexes", "Mean expression"],
        what: "Scores ligand–receptor pairs against the curated CellPhoneDB v5 database, which models multi-subunit complexes, and tests them by permuting cell-type labels.",
        when: "A solid, widely-used default for human data when you want a curated reference and proper handling of receptor/ligand complexes.",
        limits: "Human-focused; permutation p-values need enough cells per population; mean-expression scoring can favour broadly expressed genes." },
      { name: "CellChat", icon: "ti-messages", logo: "/static/img/tools/cellchat.png", tagline: "Permutation test and signalling-pathway annotation.", features: ["Complexes", "p-values", "Pathways"],
        what: "Estimates communication probability using CellChatDB v2 (ligand/receptor co-factors), aggregates pairs into signalling pathways, and assigns permutation p-values.",
        when: "When you want pathway-level summaries and statistical significance, and your populations have enough cells.",
        limits: "Slower; drops populations with fewer than ~10 cells; results reflect the CellChatDB annotations." },
      { name: "NATMI", icon: "ti-chart-dots-3", tagline: "Expression weights and specificity of each pair; fast.", features: ["Specificity", "Single genes"],
        what: "Weights each ligand–receptor pair by expression and by specificity (how restricted it is to a sender–receiver pair), at the gene level, using connectomeDB2020.",
        when: "For fast, specificity-based ranking of gene-level pairs when multi-subunit complexes are not essential.",
        limits: "No complex modelling; no significance test; gene-level pairs only." },
      { name: "SingleCellSignalR", icon: "ti-antenna-bars-5", logo: "/static/img/tools/singlecellsignalr.png", tagline: "Simple regularised ligand–receptor score.", features: ["LRscore", "Single genes"],
        what: "Computes a simple regularised ligand–receptor score (LRscore) from mean expression against the LRdb database, at the gene level.",
        when: "As a quick, interpretable complementary score alongside other methods.",
        limits: "No complexes; no permutation test; gene-level pairs (SingleCellSignalR v1)." },
      { name: "ICELLNET", icon: "ti-shield-check", tagline: "Strictly curated, human-validated interactions; immune focus.", features: ["Complexes", "Curated families"],
        what: "Scores communication against the strictly curated, experimentally-supported ICELLNET human database, grouping interactions into families (cytokines, chemokines, checkpoints…) with complex handling.",
        when: "For human and especially immune-focused datasets where high-confidence, curated families matter.",
        limits: "Human-only; a smaller curated database means lower recall; immune-centric coverage." },
    ],
    FLOW: [
      { icon: "ti-upload", title: "Add a dataset", text: "Your .h5ad file, checked and prepared" },
      { icon: "ti-route", title: "Choose a tool", text: "One is recommended for your data" },
      { icon: "ti-player-play", title: "Run", text: "Ligand–receptor interactions between populations" },
      { icon: "ti-topology-star-3", title: "Explore", text: "Interactions, network, receptors", href: "/interactions", link: "Open interactions" },
      { icon: "ti-book", title: "Document", text: "Biology, literature, your data, receiver response", href: "/interactions", link: "Open an interaction" },
      { icon: "ti-user-check", title: "Assess", text: "Your judgement, kept with each interaction" },
    ],
    roles: [
      { key: "cell_type_column", label: "Cell types", required: true, help: "The populations used as senders and receivers." },
      { key: "sample_column", label: "Samples (optional)", required: false, help: "Donor, patient or sample. Needed, with groups, to compare conditions." },
      { key: "group_column", label: "Groups (optional)", required: false, help: "Condition, for example disease and control." },
    ],

    async init() {
      this.initPanel();
      this.loadList();
      this.loadFlow();
      getJSON("/api/datasets/files").then(d => { this.files = d; }).catch(() => {});
      await this.checkEnv();
      await this.loadMethodStatus(); this.loadNcbi();
      this.loadRecommendation();
      if (!this.env.ok) this.setupOpen = true;          // open Setup only when the Python environment is missing
      this.resumeRunning();
    },

    async loadFlow() { try { this.flow = await getJSON("/api/workflow"); } catch (e) {} },
    setupOk() { return !!(this.env.ok && this.methods && this.methods.any_ok); },
    readyText() {
      if (!this.methods || !this.methods.tools) return "";
      const t = this.methods.tools, ready = this.toolList.filter(k => t[k].ok), missing = this.toolList.filter(k => !t[k].ok);
      return ready.map(k => t[k].message).join(", ") + (missing.length ? " (not installed: " + missing.join(", ") + ")" : "");
    },
    activeHasExpression() { const d = this.activeDataset(); return !!(d && d.has_expression); },
    exampleDataset() { return (this.list || []).find(d => d.id === "melanoma_tme") || null; },
    ownDatasets() { return (this.list || []).filter(d => d.id !== "melanoma_tme"); },
    toolReady(t) { return !!(this.methods && this.methods.tools && this.methods.tools[t] && this.methods.tools[t].ok); },
    isChosen(t) { return this.runMode === "single" ? this.singleTool === t : this.runTools.includes(t); },

    methodState(t) {
      if (this.checking[t]) return "checking";
      const m = this.methods && this.methods.tools && this.methods.tools[t];
      if (!m || m.checked === false) return "unchecked";
      return m.ok ? "ok" : "bad";
    },

    async loadMethodStatus() {
      try { this.methods = await getJSON("/api/methods/status"); this.rDraft = this.methods.rscript || this.rDraft; }
      catch (e) { this.methods = { tools: {}, any_ok: false }; }
    },

    // Load (check) a single tool on demand. The first time, this runs the slow
    // R-package load in the background; afterwards it is cached and instant.
    async checkTool(t, select) {
      if (this.checking[t]) return;
      this.checking = { ...this.checking, [t]: true };
      try {
        const r = await getJSON("/api/methods/check/" + encodeURIComponent(t));
        const tools = { ...(this.methods.tools || {}), [t]: { ok: r.ok, checked: true, message: r.message } };
        this.methods = { ...this.methods, tools, rscript: r.rscript || this.methods.rscript, any_ok: Object.values(tools).some(x => x.ok) };
        if (select && r.ok) {
          this._userChose = true;
          if (this.runMode === "single") this.singleTool = t;
          else if (!this.runTools.includes(t)) this.runTools = [...this.runTools, t];
        }
      } catch (e) {
        const tools = { ...(this.methods.tools || {}), [t]: { ok: false, checked: true, message: e.message } };
        this.methods = { ...this.methods, tools };
      }
      this.checking = { ...this.checking, [t]: false };
    },

    chooseTool(t) {
      const st = this.methodState(t);
      if (st === "checking") return;
      if (st === "ok") {
        this._userChose = true;
        if (this.runMode === "single") this.singleTool = t;
        else this.runTools = this.runTools.includes(t) ? this.runTools.filter(x => x !== t) : [...this.runTools, t];
        return;
      }
      this.checkTool(t, true);        // unchecked or not-yet-available: load it, then select if it works
    },
    flowState(i) {
      const f = this.flow || {};
      const done = [f.prepared || f.has_results, f.has_results, f.has_results, f.has_results, f.has_results, f.assessed > 0];
      if (i === 2 && f.running && f.running.length) return "busy";
      if (done[i]) return "done";
      const firstOpen = done.findIndex(x => !x);
      return i === firstOpen ? "ready" : "todo";
    },
    async loadRecommendation() {
      const d = this.activeDataset();
      if (!d) return;
      try {
        const a = await getJSON(`/api/datasets/${encodeURIComponent(d.id)}/advice`);
        this.rec = a.recommended;
        const t = a.recommended.tool;
        if (!this._userChose) {
          if (this.runMode === "single" && !this.singleTool) this.singleTool = t;   // pre-select the recommended one
          if (this.methodState(t) === "unchecked") this.checkTool(t);                // and load it in the background
        }
      } catch (e) {}
    },
    async runActive() {
      const d = this.activeDataset();
      if (d) await this.runMethods(d);
    },
    async deleteDataset(d) {
      if (!confirm(`Remove “${d.label}” from this computer?\n\nThis deletes its prepared data and results. Expert assessments already recorded stay in the database. This cannot be undone.`)) return;
      try {
        await getJSON(`/api/datasets/${encodeURIComponent(d.id)}/delete`, { method: "POST" });
        await this.loadList(); this.loadFlow();
      } catch (e) { alert(e.message); }
    },
    async runAnalysisFor(id) {
      try { await getJSON("/api/datasets/active", { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ id }) }); }
      catch (e) { /* navigate anyway; the run page shows a prompt if nothing loaded */ }
      location.href = "/run";
    },
    async openExample() {
      if (this.datasetsInfo.active !== "melanoma_tme") await this.switchDataset("melanoma_tme");
      else location.href = "/interactions";
    },
    async useExampleForRun() {
      if (this.datasetsInfo.active !== "melanoma_tme") {
        await getJSON("/api/datasets/active", { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ id: "melanoma_tme" }) });
        await this.loadList(); this.loadFlow(); this.loadRecommendation();
      }
      const el = document.getElementById("tools"); if (el) el.scrollIntoView({ behavior: "smooth" });
    },

    async toggleAdvice(id) {
      if (this.advice[id]) { const a = { ...this.advice }; delete a[id]; this.advice = a; return; }
      try {
        const a = await getJSON(`/api/datasets/${encodeURIComponent(id)}/advice`);
        this.advice = { ...this.advice, [id]: a };
      } catch (e) {}
    },
    adviceFacts(a) {
      const f = a.facts, parts = [];
      if (f.cells) parts.push(this.fmtInt(f.cells) + " cells");
      if (f.populations) parts.push(f.populations + " populations (smallest: " + f.smallest_population + " cells)");
      if (f.raw_counts !== null && f.raw_counts !== undefined) parts.push(f.raw_counts ? "raw counts" : "normalised values");
      return parts.length ? "Based on: " + parts.join(", ") + "." : "";
    },

    toolsToRun() {
      if (this.runMode === "single") {
        return this.singleTool && this.methods && this.methods.tools && this.methods.tools[this.singleTool] && this.methods.tools[this.singleTool].ok ? [this.singleTool] : [];
      }
      return this.chosenReady();
    },
    async viewRun(d, run) {
      try {
        await getJSON(`/api/datasets/${encodeURIComponent(d.id)}/view`, { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ run }) });
        if (d.id !== this.datasetsInfo.active) await this.switchDataset(d.id);
        else { this.loadList(); this.loadFlow(); }
      } catch (e) { alert(e.message); }
    },

    chosenReady() {
      if (!this.methods || !this.methods.tools) return [];
      return this.runTools.filter(t => this.methods.tools[t] && this.methods.tools[t].ok);
    },
    onPick(ev, kind) { const f = ev.target.files && ev.target.files[0]; if (f) this.handleUpload(f, kind); ev.target.value = ""; },
    onDrop(ev, kind) {
      if (kind === "h5ad") this.dragOver = false; else this.dragOverMeta = false;
      const f = ev.dataTransfer.files && ev.dataTransfer.files[0]; if (f) this.handleUpload(f, kind);
    },
    uploadFile(file, kind) {
      return new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.open("POST", `/api/datasets/upload?kind=${kind}&name=` + encodeURIComponent(file.name));
        xhr.upload.onprogress = (e) => { if (e.lengthComputable) this.upPct = Math.round((e.loaded / e.total) * 100); };
        xhr.onload = () => {
          if (xhr.status >= 200 && xhr.status < 300) { try { resolve(JSON.parse(xhr.responseText)); } catch (e) { reject(new Error("Bad server response")); } }
          else { let m = "Upload failed"; try { m = JSON.parse(xhr.responseText).detail || m; } catch (e) {} reject(new Error(m)); }
        };
        xhr.onerror = () => reject(new Error("Upload failed (network)"));
        xhr.send(file);
      });
    },
    async handleUpload(file, kind) {
      if (kind === "h5ad" && !file.name.toLowerCase().endsWith(".h5ad")) { this.msg = { text: "Please choose a .h5ad file.", error: true }; return; }
      this.busy = true; this.uploading = true; this.upPct = 0; this.msg = { text: `Uploading ${file.name}…`, error: false };
      try {
        const d = await this.uploadFile(file, kind);
        if (kind === "metadata") { this.form.metadata = d.path; this.metaName = d.name; }
        else { this.form.h5ad = d.path; this.h5adName = d.name; this.h5adSize = d.size_mb; }
        this.msg = { text: "", error: false };
      } catch (e) { this.msg = { text: e.message, error: true }; }
      this.busy = false; this.uploading = false; this.upPct = 0;
    },
    pickFile(f) {
      this.form.h5ad = f.path;
      const ann = (this.files.annotation || []).find(a => a.folder === f.folder);
      if (ann && !this.form.metadata) this.form.metadata = ann.path;
    },

    async checkMethods() {        // explicit "check all tools" (Setup); slow, so it shows a busy state
      this.checkingAll = true;
      try {
        this.methods = await getJSON("/api/methods/check"); this.rDraft = this.methods.rscript || this.rDraft;
      } catch (e) {}
      this.checkingAll = false;
    },

    async saveR() {
      try {
        this.methods = await getJSON("/api/settings/rscript", { method: "PUT", headers: { "Content-Type": "application/json" },
                                                               body: JSON.stringify({ rscript: this.rDraft }) });
      } catch (e) { alert(e.message); }
    },

    async saveToolsPython() {     // CellPhoneDB Python (e.g. wsl:/.../bin/python); blank clears it
      this.cpdbSaving = true;
      try {
        const r = await getJSON("/api/settings/tools-python", { method: "PUT", headers: { "Content-Type": "application/json" },
                                                                body: JSON.stringify({ tools_python: this.cpdbDraft }) });
        if (r.methods) this.methods = r.methods;
        if (r.ok === false) alert(r.message || "Could not use that Python.");
      } catch (e) { alert(e.message); }
      this.cpdbSaving = false;
    },

    async loadNcbi() {
      try {
        const s = await getJSON("/api/settings");
        this.ncbiEmail = s.ncbi_email || ""; this.ncbiKeySet = !!s.ncbi_api_key_set;
        this.ephemeral = s.ephemeral_mode !== false;
      } catch (e) {}
    },

    async saveEphemeral() {
      try {
        await getJSON("/api/settings/ephemeral", { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ ephemeral_mode: this.ephemeral }) });
        this.ephemeralSaved = this.ephemeral ? "Saved. Everything will be wiped at each launch." : "Saved. Your datasets and assessments will be kept.";
        setTimeout(() => { this.ephemeralSaved = ""; }, 3500);
      } catch (e) { this.ephemeralSaved = e.message; }
    },

    async saveNcbi() {
      try {
        const body = { ncbi_email: this.ncbiEmail };
        if (this.ncbiKey) body.ncbi_api_key = this.ncbiKey;
        const r = await getJSON("/api/settings/ncbi", { method: "PUT", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
        this.ncbiKeySet = !!r.ncbi_api_key_set; this.ncbiKey = ""; this.ncbiSaved = "Saved.";
        setTimeout(() => { this.ncbiSaved = ""; }, 2500);
      } catch (e) { this.ncbiSaved = e.message; }
    },

    async runMethods(d) {
      const again = d.has_results ? " Results of other tools stay available; you can switch between them." : "";
      const tools = this.toolsToRun();
      if (!confirm(`Run ${tools.join(", ")} for “${d.label}”? This can take from minutes to hours, depending on the size of the dataset.${again}`)) return;
      try {
        this.runJob = await getJSON(`/api/datasets/${encodeURIComponent(d.id)}/run`, { method: "POST",
          headers: { "Content-Type": "application/json" }, body: JSON.stringify({ tools }) });
        this.followRun();
      } catch (e) { alert(e.message); }
    },

    async resumeRunning() {
      try {
        const running = await getJSON("/api/jobs");
        const r = running.find(j => j.kind === "methods");
        if (r) { this.runJob = r; this.followRun(); }
      } catch (e) {}
    },

    async followRun() {
      while (this.runJob && this.runJob.state === "running") {
        await new Promise(r => setTimeout(r, 2000));
        try { this.runJob = await getJSON(`/api/jobs/${this.runJob.id}`); } catch (e) { break; }
      }
      this.loadList();
      this.loadFlow();
    },

    stepPct(step) {
      const m = /^(\d+)\/(\d+)/.exec(step || "");
      return m ? Math.round((parseInt(m[1], 10) - 1) / parseInt(m[2], 10) * 100) : 0;
    },

    async loadList() {
      try { const d = await getJSON("/api/datasets"); this.datasetsInfo = d; this.list = d.items; } catch (e) {}
    },

    async checkEnv() {
      try {
        const s = await getJSON("/api/settings");
        this.envDraft = s.analysis_python || "";
        this.cpdbDraft = s.tools_python || "";
        const r = await getJSON("/api/settings/test", { method: "POST" });
        this.env = { checked: true, ...r };
        if (!r.ok) this.envEdit = true;
      } catch (e) { this.env = { checked: true, ok: false, message: e.message }; }
    },

    async saveEnv() {
      this.envSaving = true;
      try {
        const r = await getJSON("/api/settings", { method: "PUT", headers: { "Content-Type": "application/json" },
                                                   body: JSON.stringify({ analysis_python: this.envDraft }) });
        this.env = { checked: true, ...r };
        if (r.ok) { this.envEdit = false; this.checkMethods(); }
      } catch (e) { this.env = { checked: true, ok: false, message: e.message }; }
      this.envSaving = false;
    },

    async toggleDetails(id) {
      if (this.details[id]) { const d = { ...this.details }; delete d[id]; this.details = d; return; }
      try {
        const [d, prov] = await Promise.all([
          getJSON(`/api/datasets/${encodeURIComponent(id)}`),
          getJSON(`/api/datasets/${encodeURIComponent(id)}/provenance`).catch(() => null),
        ]);
        this.details = { ...this.details, [id]: { ...d, provenance: prov } };
      } catch (e) {}
    },

    async doInspect() {
      this.busy = true; this.msg = { text: "Reading the file in the analysis environment…", error: false };
      try {
        this.info = await getJSON("/api/datasets/inspect", { method: "POST", headers: { "Content-Type": "application/json" },
                                                              body: JSON.stringify({ h5ad: this.form.h5ad, metadata: this.form.metadata }) });
        for (const role of ["cell_type", "sample", "group"]) {
          const c = (this.info.columns || []).find(x => x.suggested === role);
          this.form[role + "_column"] = c ? c.name : "";
        }
        if (!this.form.label) {
          const base = this.form.h5ad.replace(/["']/g, "").split(/[\\/]/).pop().replace(/\.h5ad$/i, "");
          this.form.label = base.replace(/[_-]+/g, " ");
          this.autoId();
        }
        this.form.matrix = (this.info.counts && this.info.counts.looks_like_raw_counts === false && this.info.has_raw) ? "raw" : "X";
        this.msg = { text: "", error: false };
        this.step = 2;
      } catch (e) { this.msg = { text: e.message, error: true }; }
      this.busy = false;
    },

    autoId() {
      if (this.idTouched) return;
      this.form.id = (this.form.label || "").toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "").slice(0, 40);
    },

    examplesOf(name) {
      const c = (this.info.columns || []).find(x => x.name === name);
      return c ? c.examples : [];
    },

    async doPrepare() {
      this.busy = true; this.msg = { text: "", error: false };
      const f = this.form;
      const body = { ...f, layer: f.matrix.startsWith("layer:") ? f.matrix.slice(6) : "", use_raw: f.matrix === "raw" };
      try {
        this.job = await getJSON("/api/datasets", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
        this.lastId = f.id;
        this.step = 3;
        this.pollJob();
      } catch (e) { this.msg = { text: e.message, error: true }; }
      this.busy = false;
    },

    async pollJob() {
      while (this.job && this.job.state === "running") {
        await new Promise(r => setTimeout(r, 1200));
        try { this.job = await getJSON(`/api/jobs/${this.job.id}`); } catch (e) { break; }
      }
      this.loadList();
    },

    copyText(text, ev) {
      navigator.clipboard.writeText(text).then(() => {
        const b = ev.target; const old = b.textContent; b.textContent = "Copied"; setTimeout(() => b.textContent = old, 1500);
      }).catch(() => {});
    },
    totalMinutes(rec) {
      const s = (rec.steps || []).reduce((a, st) => a + (st.seconds || 0), 0);
      return s ? Math.max(1, Math.round(s / 60)) : null;
    },
    swList(rec) {
      const py = (rec.software || {}).python || {}, r = (rec.software || {}).r || {};
      return [["CellPhoneDB", py.cellphonedb], ["CellChat", r.CellChat], ["NMF", r.NMF], ["presto", r.presto],
              ["R", r.R ? r.R + (r.where === "WSL" ? " (WSL)" : "") : null], ["Python", py.python],
              ["anndata", py.anndata], ["pandas", py.pandas]].filter(x => x[1]);
    },

    resetForm2() { this.form = { ...EMPTY }; this.info = {}; this.idTouched = false; this.job = {}; this.step = 1;
      this.h5adName = ""; this.h5adSize = null; this.metaName = ""; this.upPct = 0; this.uploading = false; },
  });
}
