// Builds ../results/SPOKE_tissue_validation.pptx. Run from this folder:
//   NODE_PATH=<dir with pptxgenjs> node build.js
// Numbers that come from the edge-set comparison are read from numbers.json (make_numbers.py).
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require("/Users/justinlin/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/71c655a9-053d-4261-b9d8-be4bd1cd9036/2f0a37bc-e44a-45e5-a766-4861f324561b/skills/pptx/scripts/apply_theme.js");
const N = JSON.parse(fs.readFileSync("numbers.json", "utf8"));
const THEME = { name: "SPOKE Tissue", headFontFace: "Cambria", bodyFontFace: "Calibri",
  colors: { dk1: "1B2A41", lt1: "FFFFFF", dk2: "0F4C5C", lt2: "E8EEF2", accent1: "0F4C5C", accent2: "E36414", accent3: "5F8D4E", accent4: "9A031E", accent5: "FB8B24", accent6: "3A6EA5", hlink: "0F4C5C", folHlink: "5F8D4E" } };
const pres = new pptxgen(); pres.layout = "LAYOUT_WIDE"; pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.author = "Justin Lin"; pres.title = "Control-sample tissue identity in SPOKE";
const C = pres.SchemeColor;
const foot = (dark) => ({ text: "Control gene counts mapped to SPOKE 2025  |  Justin Lin", options: { x: 0.5, y: 7.05, w: 9, h: 0.3, fontSize: 9, color: dark ? C.background2 : C.text2, isTextBox: true, margin: 0 } });
pres.defineSlideMaster({ title: "TITLE", background: { color: THEME.colors.dk2 }, objects: [
  { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.3, w: 11.7, h: 1.6, fontSize: 40, bold: true, color: C.background1, align: "left", valign: "bottom" }, text: "Title" } },
  { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.0, w: 11.7, h: 1.6, fontSize: 18, color: C.background2, align: "left", valign: "top" }, text: "Subtitle" } } ] });
pres.defineSlideMaster({ title: "SECTION", background: { color: THEME.colors.dk2 }, objects: [
  { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.8, w: 11.7, h: 1.4, fontSize: 34, bold: true, color: C.background1, align: "left" }, text: "Section" } },
  { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.2, w: 11.7, h: 1.2, fontSize: 16, color: C.background2 }, text: "" } } ] });
pres.defineSlideMaster({ title: "CONTENT", background: { color: THEME.colors.lt1 }, objects: [
  { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.35, w: 12.3, h: 0.9, fontSize: 28, bold: true, color: C.text2, align: "left", valign: "middle", margin: 0 }, text: "Title" } },
  { text: foot(false) }, { text: { text: "", options: { x: 0, y: 0, w: 0, h: 0 } } } ], slideNumber: { x: 12.4, y: 7.05, w: 0.5, h: 0.3, fontSize: 9, color: C.text2 } });
const bullets = (items, opts = {}) => items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1, paraSpaceAfter: 6, color: C.text1, fontSize: opts.fontSize || 16 } }));
const tableOpts = (fs_) => ({ fontSize: fs_ || 11, fontFace: "Calibri", color: THEME.colors.dk1, border: { type: "solid", pt: 0.5, color: "C9D3DB" }, valign: "middle", autoPage: false });
const hdr = (cells) => cells.map(c => ({ text: c, options: { bold: true, fill: { color: THEME.colors.lt2 }, color: THEME.colors.dk2 } }));
const IMG = { liver_qq: "../datasets/OSD-245__2025upAuG/corrected_ontology/A_qq_plot.png", builds: "../results/pilot_stage2_builds_2025only.png",
  overview: "../results/tissue_overview_2025upAuG60.png", liver_full: "../datasets/OSD-245__2025/corrected_linked_2025/A_linked_nodes.png",
  liver_up: "../datasets/OSD-245__2025upAuG/corrected_linked_2025upAuG60/A_linked_nodes.png" };
for (const p of Object.values(IMG)) if (!fs.existsSync(p)) throw new Error("missing figure " + p);

// 1 title
let s = pres.addSlide({ masterName: "TITLE" });
s.addText("Does a control transcriptome point at its own tissue in SPOKE?", { placeholder: "title" });
s.addText("A validated configuration for projecting NASA OSDR mouse gene counts onto the SPOKE knowledge graph\nJustin Lin  ·  104 control sets, 31 tissues, 916 samples  ·  SPOKE 2025-03", { placeholder: "body" });

// 2 why
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Why control samples, and why anatomy", { placeholder: "title" });
s.addText(bullets([
  "Goal: use SPOKE embeddings (PSEVs) of OSDR mouse RNA-seq for downstream models of spaceflight effects.",
  "Before trusting any flight result, the projection must recover something we already know: the tissue the sample came from.",
  "Control samples carry no spaceflight signal, so their Anatomy-node ranking is a built-in positive control.",
  "Every dataset is judged on its own: no dataset is compared with another.",
  "Scope: every mouse bulk RNA-seq study on OSDR with GeneLab counts and a control group: 94 studies, 104 control sets, 31 tissues, 916 samples (83 spaceflight ground controls, 21 ground-analog baselines)."]), { x: 0.5, y: 1.4, w: 7.6, h: 5.4, valign: "top", isTextBox: true });
s.addShape(pres.ShapeType.roundRect, { x: 8.5, y: 1.6, w: 4.3, h: 4.6, fill: { color: THEME.colors.lt2 }, line: { color: THEME.colors.lt2 }, rectRadius: 0.15 });
s.addText([{ text: "104", options: { fontSize: 54, bold: true, color: C.accent1, breakLine: true } }, { text: "control datasets", options: { fontSize: 14, color: C.text2, breakLine: true } },
  { text: "31", options: { fontSize: 54, bold: true, color: C.accent1, breakLine: true } }, { text: "tissues", options: { fontSize: 14, color: C.text2, breakLine: true } },
  { text: "13,895", options: { fontSize: 30, bold: true, color: C.accent1, breakLine: true } }, { text: "Anatomy nodes in SPOKE, 60 with enough gene edges to test", options: { fontSize: 12, color: C.text2 } }],
  { x: 8.7, y: 1.8, w: 3.9, h: 4.2, align: "center", valign: "middle", isTextBox: true, margin: 0 });

// 3 pipeline
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Pipeline: counts to a tissue verdict", { placeholder: "title" });
const steps = [["1", "Counts → weights", "GeneLab RSEM counts, control columns only. CPM per sample; keep genes with CPM ≥ 1; weight = log2(CPM + 1)."],
  ["2", "Mouse → SPOKE genes", "Ensembl → mouse Entrez (GeneLab annotation) → human Entrez (HomoloGene 68) → SPOKE Gene node."],
  ["3", "PSEV projection", "Per gene: z-score its PSEV over the dataset's genes, rank the 547,869 nodes within the gene. Per sample: weighted mean of ranks over all genes."],
  ["4", "Baseline correction", "Re-run with every gene weighted 1 (the uniform profile); subtract it scaled by the sample's mean weight."],
  ["5", "Rank & test", "Primary: rank of the tissue's own node among the 60 gene-linked Anatomy nodes. Secondary: UBERON family vs all Anatomy nodes (QQ plot, Mann-Whitney U, KS). Fallback: the tissue's cell type among 2,745 CellType nodes."]];
steps.forEach((st, i) => { const y = 1.45 + i * 1.08;
  s.addShape(pres.ShapeType.ellipse, { x: 0.6, y: y + 0.08, w: 0.6, h: 0.6, fill: { color: THEME.colors.accent1 }, line: { color: THEME.colors.accent1 } });
  s.addText(st[0], { x: 0.6, y: y + 0.08, w: 0.6, h: 0.6, align: "center", valign: "middle", fontSize: 16, bold: true, color: C.background1, isTextBox: true, margin: 0 });
  s.addText([{ text: st[1], options: { bold: true, fontSize: 15, color: C.text2, breakLine: true } }, { text: st[2], options: { fontSize: 12.5, color: C.text1 } }], { x: 1.45, y: y, w: 11.2, h: 0.95, valign: "middle", isTextBox: true, margin: 0 }); });

// 4 section
s = pres.addSlide({ masterName: "SECTION" }); s.addText("Three findings that shaped the configuration", { placeholder: "title" });
s.addText("1. The raw ranking is mostly baseline   ·   2. Few Anatomy nodes carry gene links   ·   3. Edge type decides everything", { placeholder: "body" });

// 5 finding 1
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Finding 1: the raw ranking is mostly the graph, not the sample", { placeholder: "title" });
s.addText(bullets([
  "A liver dataset's raw Anatomy ranking correlates at Spearman 0.93 with the ranking you get when every gene is weighted equally.",
  "Liver, kidney, spleen, retina, heart, cerebellum and muscle datasets agree with one another at a median of 0.87; one shared component carries 71 % of their variance. The same five kidney-cortex nodes head every one of them.",
  "Cause: all weights are positive, so a node's score is dominated by its average rank across genes, a fixed property of the graph.",
  "Fix (within each dataset): subtract the score the node gets when every gene is weighted equally. What remains is the covariance between expression and the node's gene links."], { fontSize: 14.5 }), { x: 0.5, y: 1.4, w: 6.6, h: 4.2, valign: "top", isTextBox: true });
s.addTable([hdr(["Liver dataset OSD-245", "raw rank", "corrected rank"]),
  ["right lobe of liver, of 13,895", "161", "9"], ["right lobe of liver, of 60 linked", "4", "1"], ["cortex of kidney, of 13,895", "1", "2,882"]],
  { x: 0.5, y: 5.65, w: 6.6, colW: [3.4, 1.6, 1.6], rowH: 0.28, ...tableOpts(11), align: "left" });
s.addImage({ path: IMG.liver_qq, x: 7.5, y: 1.35, w: 5.3, h: 5.47 });
s.addText("Liver family (70 UBERON nodes) after correction: median percentile 0.3, all 19 samples", { x: 7.5, y: 6.78, w: 5.3, h: 0.25, fontSize: 9.5, color: C.text2, isTextBox: true, margin: 0, align: "center" });

// 6 finding 2
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Finding 2: only a few dozen Anatomy nodes carry gene links", { placeholder: "title" });
s.addText(bullets([
  "SPOKE loads all of UBERON (13,895 Anatomy nodes), but Anatomy–Gene edges come from human expression atlases with small tissue panels: 92 nodes have any such edge, 72 have up-regulation edges, 60 have at least 200 of them.",
  "The panel is GTEx-style ('right lobe of liver', 'upper lobe of left lung', 'skin of leg'). Retina, thymus, bone marrow, eye and placenta have no linked node at all.",
  "Every other node scores second-hand through the ontology. The honest test is: among the linked nodes, where does the tissue's own node rank?",
  `On the full graph that test exposes a connectivity artifact: a node's rank in its own tissue is predicted by how many edges it has (Spearman ${N.spearman_full}). The most lightly linked nodes top almost every dataset; liver ranks ${N.liver_full} of ${N.n_linked_full}, spleen ${N.spleen_full}, blood ${N.blood_full}.`], { fontSize: 15 }), { x: 0.5, y: 1.4, w: 12.3, h: 5.3, valign: "top", isTextBox: true });

// 7 finding 3
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Finding 3: edge type is the cause, and 'upregulates' edges fix it", { placeholder: "title" });
s.addImage({ path: IMG.builds, x: 0.5, y: 1.35, w: 6.6, h: 5.49 });
s.addText(bullets([
  "The PSEV walk uses a binary adjacency: an 'expresses', 'downregulates' or 'upregulates' edge counts the same.",
  "1.58 M blanket 'expresses' edges (blood → 18,879 of 19,072 genes), 494 k 'downregulates', 86 k 'upregulates' (log2 fold-change ≥ 1, median 1.6).",
  `Full graph: lightly linked nodes win everywhere (connectivity Spearman ${N.spearman_full}).`,
  "Without 'expresses' but with 'downregulates': fails everywhere.",
  "Housekeeping genes removed, full graph: no change. The problem is in the edges, not the genes.",
  `'Upregulates' only: every testable tissue's own node in the top three, connectivity correlation gone (${N.spearman_up}), no node wins everywhere.`], { fontSize: 13.5 }), { x: 7.4, y: 1.4, w: 5.5, h: 5.4, valign: "top", isTextBox: true });

// 8 section + 9 configuration
s = pres.addSlide({ masterName: "SECTION" }); s.addText("The validated configuration", { placeholder: "title" }); s.addText("Same pipeline, different edge set", { placeholder: "body" });
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Validated configuration", { placeholder: "title" });
const cfg = [["Graph", "SPOKE 2025-03 snapshot, gene PSEVs built with only the 85,530 UPREGULATES Anatomy–Gene edges; every other edge type unchanged (547,869 nodes)"],
  ["Weights", "log2(CPM + 1) per gene and sample, CPM ≥ 1 filter, mouse → human by HomoloGene build 68"],
  ["Score", "rank-based PSEV projection with within-dataset baseline correction (uniform-weight profile subtracted)"],
  ["Anatomy test", "rank of the tissue's own node among the 60 Anatomy nodes with ≥ 200 up-regulation edges; CLEAR = 1st, PARTIAL = top 7"],
  ["Cell-type route", "for tissues without a linked Anatomy node: rank of the matching CellType node among 2,745 (59 carry gene edges)"],
  ["Family test", "UBERON-defined tissue family vs all other Anatomy nodes: QQ plot, one-sided Mann-Whitney U and KS, per-sample consistency"],
  ["What was rejected", "z-score projection, hub-matched permutation null, housekeeping-gene removal: none changed the failing tissues"]];
s.addTable([hdr(["Component", "Choice"]), ...cfg.map(r => [{ text: r[0], options: { bold: true, color: THEME.colors.dk2 } }, r[1]])],
  { x: 0.5, y: 1.4, w: 12.3, colW: [2.4, 9.9], rowH: 0.62, ...tableOpts(12.5), align: "left" });

// 10 results
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Results: 104 datasets on the validated configuration", { placeholder: "title" });
s.addImage({ path: IMG.overview, x: 0.5, y: 1.3, w: 4.1, h: 5.57 });
s.addTable([hdr(["Own-node outcome", "datasets"]), ["ranked 1st (CLEAR)", "50"], ["2nd–7th (PARTIAL)", "32"], ["below 7th (NONE)", "2"], ["no linked Anatomy node", "20"]],
  { x: 5.0, y: 1.4, w: 4.6, colW: [3.2, 1.4], rowH: 0.34, ...tableOpts(12), align: "left" });
s.addText(bullets([
  "1st in every dataset: liver (14/14), all five hindlimb muscles (20/20), cerebrum, brain, white adipose, blood; kidney 4/6 (the other two 2nd); cerebellum 2/3; skin 1/7 + 6 at 2nd–4th.",
  "Consistent near-misses behind a related node: spleen 2nd (blood) in 6/8, heart 2nd (gastrocnemius) in 6/6, adrenal 2nd–3rd (aorta), lung 4th–6th (aorta, arteries), hippocampus 4th (prefrontal cortex).",
  "Failures: one colon set and brown adipose.",
  "Ground-analog baselines behave like spaceflight ground controls (11 CLEAR / 6 PARTIAL / 4 no node)."], { fontSize: 12.5 }), { x: 5.0, y: 3.4, w: 7.9, h: 3.4, valign: "top", isTextBox: true });

// 11 three groups
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Tissue coverage in three groups", { placeholder: "title" });
s.addTable([hdr(["Group", "Tissues", "How validated"]),
  [{ text: "Validated directly", options: { bold: true, color: THEME.colors.accent3 } }, "liver, kidney, heart, spleen, lung, adrenal gland, skin, brain, cerebrum, cerebellum, hippocampus, white adipose, blood, mammary gland, gastrocnemius, soleus, quadriceps, tibialis anterior, extensor digitorum longus, colon (weak)", "own Anatomy node among the 60 linked nodes"],
  [{ text: "Validated by cell type", options: { bold: true, color: THEME.colors.accent5 } }, "retina (photoreceptor cell, rank 1 in 4/5, 2 in the fifth), bone marrow (monocyte, 1), placenta (syncytiotrophoblast, 1), cecum (Paneth cell, 2); thymus weak (T cell, median 17 of 2,745); eye marginal (rod cell, 15–23)", "own CellType node among 2,745"],
  [{ text: "Not testable", options: { bold: true, color: THEME.colors.accent4 } }, "bone, mandible, temporal bone, optic nerve", "no Anatomy or CellType node with gene edges; only ontology-neighbour proxies"]],
  { x: 0.5, y: 1.4, w: 12.3, colW: [2.3, 7.0, 3.0], rowH: 1.05, ...tableOpts(12.5), align: "left" });
s.addText("27 of 31 tissues have a direct within-dataset test; thymus, eye, colon and brown adipose carry an explicit 'weak' label.", { x: 0.5, y: 5.2, w: 12.3, h: 0.5, fontSize: 13, color: C.text2, isTextBox: true, margin: 0 });
s.addText("Why the gaps: expression edges come from human atlases (GTEx-style sites); GTEx never sampled thymus, retina, bone marrow or placenta, and no atlas edge exists for bone.", { x: 0.5, y: 5.75, w: 12.3, h: 0.7, fontSize: 12, color: C.text2, isTextBox: true, margin: 0 });

// 12 liver example
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Example: liver (OSD-245, 19 ground controls)", { placeholder: "title" });
s.addImage({ path: IMG.liver_full, x: 0.5, y: 1.55, w: 5.0, h: 5.2 });
s.addImage({ path: IMG.liver_up, x: 7.6, y: 1.55, w: 5.2, h: 5.2 });
s.addText(`Full graph, all Anatomy–Gene edges: best liver node ${N.liver_full_ordinal} of ${N.n_linked_full}`, { x: 0.5, y: 1.3, w: 5.0, h: 0.3, fontSize: 12, bold: true, color: C.accent4, isTextBox: true, margin: 0, align: "center" });
s.addText("Upregulates edges only: right lobe of liver 1st of 60", { x: 7.6, y: 1.3, w: 5.2, h: 0.3, fontSize: 12, bold: true, color: C.accent3, isTextBox: true, margin: 0, align: "center" });
s.addText([{ text: "Same samples,\nsame weights,\nsame scoring.", options: { fontSize: 16, bold: true, color: C.text2, breakLine: true } }, { text: "\nOnly the edge set differs.", options: { fontSize: 13, color: C.text1 } }], { x: 5.6, y: 3.0, w: 1.9, h: 2.4, align: "center", valign: "middle", isTextBox: true, margin: 0 });

// 13 edge-set table
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("The same ten tissues on four edge sets of the 2025 graph", { placeholder: "title" });
s.addTable([hdr(["tissue", `all edges (${N.n_linked_full} linked)`, `no 'expresses' (${N.n_linked_noAeG})`, `housekeeping removed (${N.n_linked_nohk})`, `'upregulates' only (${N.n_linked_up})`]),
  ...N.stage2_rows.map(r => r.map((v, j) => (j === 4 && (v === "1" || v === "2" || v === "3")) ? { text: v, options: { bold: true, color: THEME.colors.accent3 } } : v)),
  [{ text: "Spearman(link weight, own rank)", options: { bold: true } }, N.spearman_full, N.spearman_noAeG, N.spearman_nohk, N.spearman_up]],
  { x: 0.5, y: 1.4, w: 12.3, colW: [2.9, 2.35, 2.35, 2.35, 2.35], rowH: 0.36, ...tableOpts(12), align: "center" });
s.addText("Rank of the tissue's own gene-linked node among that edge set's linked nodes (pooled over the dataset's samples; 'no node' = the tissue has no linked node on that edge set). Only the Anatomy–Gene edge types change between columns.", { x: 0.5, y: 5.75, w: 12.3, h: 0.8, fontSize: 12, color: C.text2, isTextBox: true, margin: 0 });

// 14 limitations
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Limitations and caveats", { placeholder: "title" });
s.addText(bullets([
  "Coverage follows the human atlases: 60 testable Anatomy nodes, 59 cell types. Mouse tissues are validated through human counterparts.",
  "Residual hubs: gastrocnemius is in the top four of 43 datasets, sinoatrial myocytes in a quarter of cell-type rankings; they edge out heart and adrenal gland.",
  "Family definitions are ontology-based (UBERON is-a / part-of) and still mix parenchyma with vessels and primordia for some organs.",
  "Small sets: several tissues have one dataset; OSD-270 heart has two samples; OSD-289 'GC' includes vivarium controls.",
  "RR-1 thymus and liver counts share a strong technical factor (replicates pair across tissues at r ≈ 0.99).",
  "Direction is lost: 'downregulates' edges had to be dropped because the walk is unsigned. A signed or two-node representation would recover them."], { fontSize: 14.5 }), { x: 0.5, y: 1.4, w: 12.3, h: 5.4, valign: "top", isTextBox: true });

// 15 next steps
s = pres.addSlide({ masterName: "CONTENT" }); s.addText("Next steps", { placeholder: "title" });
const nx = [["Freeze", "Configuration, tissue table, pipeline and methods are committed; the standalone repository runs the check on any counts table."],
  ["Benchmark", "Tissue classifier on corrected SPOKE vectors vs gene expression, held out by study: tests cross-study transfer and whether tissue identity survives outside the Anatomy nodes."],
  ["Flight", "Project the flight samples with the same pipeline, confirm the tissue still ranks first, then test flight vs ground within study and tissue, starting from known effects."],
  ["Graph", "Background experiment: rebuild with only up-regulation edges of log2 fold-change ≥ 2; add mouse atlas edges only if thymus or bone become important."]];
nx.forEach((n, i) => { const x = 0.5 + i * 3.1;
  s.addShape(pres.ShapeType.roundRect, { x, y: 1.6, w: 2.9, h: 4.4, fill: { color: THEME.colors.lt2 }, line: { color: THEME.colors.lt2 }, rectRadius: 0.12 });
  s.addText([{ text: String(i + 1), options: { fontSize: 36, bold: true, color: C.accent1, breakLine: true } }, { text: n[0], options: { fontSize: 18, bold: true, color: C.text2, breakLine: true } }, { text: "\n" + n[1], options: { fontSize: 12.5, color: C.text1 } }], { x: x + 0.2, y: 1.8, w: 2.5, h: 4.0, valign: "top", isTextBox: true, margin: 0 }); });
s.addText("Methods and every number: GC Counts Check/README.md and METHODS.md  ·  pipeline: github.com/justinytlin/SPOKE-GeneCounts-Anatomy-Pipeline  ·  project: github.com/justinytlin/NASA-OSDR-SPOKE-Integration", { x: 0.5, y: 6.25, w: 12.3, h: 0.4, fontSize: 11, color: C.text2, isTextBox: true, margin: 0 });

(async () => { await pres.writeFile({ fileName: "../results/SPOKE_tissue_validation.pptx" }); await applyTheme("../results/SPOKE_tissue_validation.pptx", THEME); console.log("deck written"); })();
