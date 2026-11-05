/* Help page: glossary and data sources. */
function helpApp() {
  const TERMS = [
    ["Interaction", "A sender population, a receiver population, a ligand and a receptor. Complexes are kept whole."],
    ["Tool", "One of the cell–cell communication methods: CellPhoneDB, NATMI, CellChat, SingleCellSignalR or ICELLNET. Each is used with its own criteria and database; you choose which to run."],
    ["Reported by", "The tools whose results include the interaction."],
    ["Raw score", "The score a tool gives to its prediction. Its meaning differs between tools, so scores are never compared between tools. “What these numbers mean” in the panel defines each one."],
    ["Rank in tool", "Where the raw score falls within that tool’s own results for this dataset (100 = highest). A rank, not a probability; ranks of different tools are not compared."],
    ["Complex", "A receptor (or ligand) made of several proteins, for example CD44+CD74. All subunits are kept and shown."],
    ["Small molecule", "A non-protein ligand, such as prostaglandin E2, estimated from the expression of the enzyme that makes it."],
    ["Sender / receiver", "The population expressing the ligand, and the population expressing the receptor."],
    ["Detected in", "Share of cells of a population with at least one count of the gene. For a complex, the share of cells in which all subunits are detected together."],
    ["Pair-specific evidence", "A publication that reports something about this ligand–receptor pair, not only about each molecule."],
    ["Co-mention only", "A publication that mentions both molecules without reporting their interaction. It is not evidence of the interaction."],
    ["Receptor signalling in general", "Curated downstream effects of the receptor (SIGNOR) or pathways that include it (Reactome), whatever ligand binds it. Not shown to be active in your data."],
    ["Receiver response", "In the Response tab: within the receiver population, cells with the receptor are compared with cells without it for pathway activity (PROGENy), transcription factor activity (CollecTRI) and target genes. It shows functional consistency, not causality."],
    ["Fold enrichment", "In Receptors: how many times more often a receptor appears in the chosen set of interactions than among all interactions considered."],
    ["FDR", "False discovery rate: the p-value of the Fisher exact test corrected for testing all receptors (Benjamini–Hochberg). Small values mean the over-representation is unlikely to be chance."],
    ["Overlap", "In Tools: interactions reported by both tools, divided by those reported by either."],
    ["Analysis record", "The file written by every run of the tools, with software and database versions, parameters, checksums and timings."],
  ];
  return Object.assign(docPanel(), {
    q: "", sources: [],
    async init() { this.initPanel(); try { this.sources = await getJSON("/api/sources"); } catch (e) {} },
    terms() { const q = this.q.trim().toLowerCase(); return q ? TERMS.filter(t => (t[0] + " " + t[1]).toLowerCase().includes(q)) : TERMS; },
  });
}
