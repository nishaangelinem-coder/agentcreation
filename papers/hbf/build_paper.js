// Builds the HBF manuscript as .docx (docx-js). Run: node build_paper.js
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  AlignmentType, WidthType, ShadingType, BorderStyle, Header, Footer, PageNumber,
  TabStopType, LevelFormat, HeadingLevel,
} = require("docx");

const FONT = "Times New Roman";
const TEXT_W = 9026; // A4 width minus 1" margins, DXA

// ---------- inline markup: **bold**, *italic*, _{sub}, ^{sup} ----------
function runs(text, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|_\{[^}]+\}|\^\{[^}]+\})/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index), ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true, ...base }));
    else if (t.startsWith("*")) out.push(new TextRun({ text: t.slice(1, -1), italics: true, ...base }));
    else if (t.startsWith("_{")) out.push(new TextRun({ text: t.slice(2, -1), subScript: true, ...base }));
    else out.push(new TextRun({ text: t.slice(2, -1), superScript: true, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), ...base }));
  return out;
}

const body = [];
const P = (t, o = {}) => body.push(new Paragraph({ children: runs(t), alignment: AlignmentType.JUSTIFIED,
  spacing: { after: 120, line: 288 }, indent: o.noIndent ? undefined : { firstLine: 360 }, ...o.p }));
const H1 = (t) => body.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun(t)] }));
const H2 = (t) => body.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun(t)] }));
const BUL = (items) => items.forEach((t) => body.push(new Paragraph({ numbering: { reference: "bul", level: 0 },
  children: runs(t), alignment: AlignmentType.JUSTIFIED, spacing: { after: 80, line: 288 } })));
const NUM = (items, ref) => items.forEach((t) => body.push(new Paragraph({ numbering: { reference: ref, level: 0 },
  children: runs(t), alignment: AlignmentType.JUSTIFIED, spacing: { after: 80, line: 288 } })));
const EQ = (t, n) => body.push(new Paragraph({
  tabStops: [{ type: TabStopType.CENTER, position: TEXT_W / 2 }, { type: TabStopType.RIGHT, position: TEXT_W }],
  spacing: { before: 80, after: 140 },
  children: [new TextRun("\t"), ...runs(t, { italics: true }), new TextRun({ text: `\t(${n})` })] }));

function pngSize(f) { const b = fs.readFileSync(f); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
let figN = 0, tabN = 0;
function FIG(file, caption, widthIn = 6.2) {
  figN++;
  const [w, h] = pngSize(file); const W = widthIn * 96;
  body.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 }, keepNext: true,
    children: [new ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: W, height: W * h / w },
      altText: { title: `Figure ${figN}`, description: caption, name: `fig${figN}` } })] }));
  body.push(new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { after: 200 },
    children: [new TextRun({ text: `Fig. ${figN}. `, bold: true, size: 20 }), ...runs(caption, { size: 20 })] }));
}
const border = { style: BorderStyle.SINGLE, size: 4, color: "808080" };
function TBL(caption, header, rows, widths) {
  tabN++;
  const total = widths.reduce((a, b) => a + b, 0);
  body.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 80 }, keepNext: true,
    children: [new TextRun({ text: `Table ${tabN}. `, bold: true, size: 20 }), ...runs(caption, { size: 20 })] }));
  const cell = (t, i, head) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    borders: { top: border, bottom: border, left: border, right: border },
    shading: head ? { fill: "D9E2EC", type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 50, bottom: 50, left: 90, right: 90 },
    children: [new Paragraph({ alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
      children: runs(String(t), { size: 18, bold: head || undefined }) })] });
  body.push(new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: widths,
    rows: [new TableRow({ tableHeader: true, children: header.map((t, i) => cell(t, i, true)) }),
      ...rows.map((r) => new TableRow({ cantSplit: true, children: r.map((t, i) => cell(t, i, false)) }))] }));
  body.push(new Paragraph({ spacing: { after: 120 }, children: [] }));
}

// =====================================================================
// FRONT MATTER
// =====================================================================
const front = [
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 },
    children: [new TextRun({ text: "High Bandwidth Flash for Large Language Model Inference: Prefill–Decode Mechanics, System Integration, and the Memory-Hierarchy Architecture Debate", bold: true, size: 32 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 },
    children: [new TextRun({ text: "M. Nisha Angeline", size: 24, bold: true }), new TextRun({ text: "*", superScript: true, size: 24 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
    children: [new TextRun({ text: "Professor and Head, Department of Electronics and Communication Engineering,", italics: true, size: 21 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
    children: [new TextRun({ text: "Velalar College of Engineering and Technology, Thindal, Erode, Tamil Nadu, India", italics: true, size: 21 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 300 },
    children: [new TextRun({ text: "*Corresponding author: nishavlsidesign@gmail.com", size: 20 })] }),
];
body.push(...front);

body.push(new Paragraph({ spacing: { after: 80 }, children: [new TextRun({ text: "Abstract", bold: true, size: 24 })] }));
P("The memory system, not arithmetic, now bounds the cost and latency of large language model (LLM) inference. High Bandwidth Memory (HBM) supplies terabytes per second but only tens of gigabytes per stack, forcing trillion-parameter mixture-of-experts (MoE) models to be sharded across many accelerator packages purely for capacity. High Bandwidth Flash (HBF) — vertically stacked 3D NAND dies joined by through-silicon vias to a logic base die and exposed through an HBM-like wide interface — has been proposed to deliver roughly an order of magnitude more capacity per stack at HBM-class read bandwidth, at the price of microsecond-scale access latency and finite write endurance. This paper develops a first-principles analysis of where HBF helps and where it hurts. We (i) derive the operational-intensity, capacity, and concurrency equations that govern the prefill and decode phases on a tiered HBM–HBF package, including the expert-touch probability that makes MoE decode traffic batch-dependent; (ii) propose a reference integration architecture covering the physical package, a read-optimised flash controller on the base die, a data-placement policy and the serving-runtime changes HBF requires; (iii) evaluate four package configurations on dense 70B and 1T-parameter MoE workloads with an open, parameterised analytical model under a 50 ms time-per-output-token service-level objective; and (iv) structure the ongoing architecture debate among HBM scale-out, on-package flash, CXL-attached pools, processing-near-memory and SRAM-centric dataflow designs. The model shows that HBF never improves per-package throughput for models that already fit in HBM, but for a 1T MoE at 32k context it enables single-package serving that HBM cannot provide and yields 1.75× (0.8 TB/s per HBF stack) to 4.2× (1.6 TB/s) higher per-package throughput at a four-package footprint. Endurance analysis shows that HBF is viable for weights and write-once prefix caches (≥8 years at one full rewrite per day) but not for the mutable KV cache (≈9 days at 1 GB/s on TLC). Sustaining HBM-class bandwidth at 20 µs array latency requires ≈1,000 concurrent 16 KiB reads per stack. We conclude that the most defensible near-term role for HBF is as a weight tier for capacity-constrained MoE deployments and for compute-bound prefill nodes in disaggregated serving, while HBM remains the home of the KV cache.", { noIndent: true });
body.push(new Paragraph({ spacing: { after: 240 }, alignment: AlignmentType.JUSTIFIED, children: [
  new TextRun({ text: "Keywords: ", bold: true }), new TextRun({ text: "High Bandwidth Flash (HBF); High Bandwidth Memory (HBM); 3D NAND; large language model inference; prefill and decode; KV cache; mixture of experts; memory hierarchy; heterogeneous integration; through-silicon via; disaggregated serving; roofline model." })] }));

// =====================================================================
H1("1. Introduction");
P("Transformer-based large language models [1], [2] have moved from research artefacts to infrastructure that serves billions of requests per day. Unlike training, which can be organised to reuse each weight across thousands of tokens, autoregressive inference generates one token per sequence per step. In this *decode* phase every generated token requires streaming the active model weights and the per-sequence key–value (KV) cache from memory, so the achievable throughput is set by memory bandwidth and the achievable batch by memory capacity [3], [4]. Gholami *et al.* [3] documented that peak accelerator FLOP/s has grown far faster than DRAM bandwidth and interconnect bandwidth over two decades; the resulting \"memory wall\" is now the dominant constraint for serving.");
P("High Bandwidth Memory (HBM) [23], [24] has been the industry's answer: stacks of DRAM dies connected by through-silicon vias (TSVs) and placed beside the accelerator on a silicon interposer provide on the order of 1 TB/s per stack. Yet HBM capacity grows slowly — a 12-high HBM3E stack holds 36 GB — and its cost per bit is several times that of commodity DRAM. Meanwhile frontier models have grown to hundreds of billions or trillions of parameters, increasingly as sparse mixture-of-experts (MoE) networks [15]–[18], [44] whose *active* parameter count per token is small but whose *total* parameter count must still be resident. A 1T-parameter MoE stored in 8-bit precision occupies ≈1 TB before any KV cache is allocated, which exceeds the 288 GB of an eight-stack HBM3E package by more than 3×. Deployments therefore shard such models across four to sixteen packages for *capacity* rather than for compute, paying for idle FLOP/s, interconnect and power.");
P("High Bandwidth Flash (HBF) attacks this capacity gap directly. Publicly announced in 2025 [29] and subsequently the subject of a standardisation memorandum between a NAND and an HBM vendor [30], HBF stacks 3D NAND core dies [26], [28] with TSVs onto a logic base die and exposes them to the host through an HBM-like wide, parallel interface. The stated goal is HBM-class *read* bandwidth with 8–16× the capacity per stack at comparable cost. NAND flash, however, differs from DRAM in three fundamental ways: (a) random read latency is tens of microseconds rather than tens of nanoseconds; (b) writes are page-granular, require prior block erasure and are orders of magnitude slower than reads; and (c) each cell tolerates only a finite number of program/erase (P/E) cycles [27]. Whether HBF is a revolution for inference or a niche capacity tier therefore depends on a careful accounting of *which* data structures are read-mostly, *which* phase of inference is bandwidth- versus compute-bound, and *how* latency can be hidden.");
P("This paper provides that accounting. Its contributions are:");
NUM([
  "**Phase-resolved analysis.** We derive closed-form expressions for the operational intensity of prefill and decode on a tiered HBM–HBF package, the minimum prefill chunk that keeps flash-resident weights compute-bound, the expert-touch probability that makes MoE weight traffic grow with batch size, and the Little's-law concurrency that a flash stack must sustain to reach its rated bandwidth (Section 3).",
  "**Reference integration architecture.** We specify a physical, controller-level and software-level integration of HBF, including a data-placement policy for weights, KV cache, prefix caches, adapters and activations, and quantify the endurance budget that rules each placement in or out (Section 4).",
  "**Quantitative evaluation.** Using an open, parameterised analytical model, we compare HBM-only and three HBM/HBF mixes for dense 70B and 1T-parameter MoE models across package counts, context lengths and HBF bandwidth projections under a 50 ms time-per-output-token (TPOT) service-level objective (Sections 5–6).",
  "**Structured architecture debate.** We contrast five competing visions of the inference memory hierarchy and identify the regimes in which each is preferable, concluding with a heterogeneous \"HBF-for-prefill, HBM-for-decode\" fleet proposal and an agenda of open research problems (Sections 7–8).",
], "contrib");

// =====================================================================
H1("2. Background and Related Work");
H2("2.1 The two phases of autoregressive inference");
P("A request to an LLM proceeds in two phases. In **prefill**, the entire prompt of *L* tokens is processed in parallel: every linear layer performs a matrix–matrix product between an *L*-row activation matrix and the weight matrix, and the keys and values of every prompt token are written into the KV cache. In **decode**, the model generates output tokens one at a time; each step multiplies a single activation vector per sequence by every active weight matrix and attends over the full KV cache of that sequence. Batching *B* concurrent sequences turns decode's matrix–vector products into thin matrix–matrix products, raising arithmetic intensity proportionally to *B*, but each sequence's attention still reads its own KV cache, which does not amortise across the batch [4], [6].");
P("The two phases therefore stress hardware differently: prefill is typically compute-bound and determines time-to-first-token (TTFT), whereas decode is memory-bound and determines time-per-output-token (TPOT). Production systems interleave them through iteration-level scheduling [7] and chunked prefill [8], or physically separate them onto different machine pools [9]–[11]. These system choices, as we show in Section 7, interact strongly with the choice of memory technology.");
H2("2.2 High Bandwidth Memory");
P("HBM stacks DRAM dies on a base logic die, connected by thousands of TSVs and a 1024-bit (HBM3) or 2048-bit (HBM4) wide interface routed through a silicon interposer or bridge [23], [24]. Its key advantage is bandwidth density at the package shoreline; its limitations are capacity per stack, cost per bit, thermal density of the stack, and the limited \"beachfront\" around the accelerator die that caps the number of stacks — typically six to eight. Fine-grained DRAM organisations [25] and processing-in-memory variants [36] have been proposed to improve energy and effective bandwidth, but neither alters the fundamental capacity-per-stack limit of DRAM cells.");
H2("2.3 NAND flash and High Bandwidth Flash");
P("3D NAND stores charge in vertically stacked cells formed in a single high-aspect-ratio etch [28]; contemporary parts exceed 200 word-line layers and store three (TLC) or four (QLC) bits per cell, yielding areal densities an order of magnitude above DRAM. Reading a page requires sensing a string through the vertical channel, which takes tens of microseconds (tR); programming takes hundreds of microseconds and erasure milliseconds. Reliability is managed by strong ECC (typically LDPC), read-retry, wear levelling and periodic refresh to counter retention loss and read disturb [27].");
P("HBF, as described in public disclosures [29], [30], re-architects NAND for bandwidth rather than for storage density alone. Its principal ingredients are: (i) TSV-connected stacking of many NAND core dies above a logic base die, analogous to HBM; (ii) partitioning of each die into a large number of independently addressable sub-arrays or planes so that thousands of reads can proceed concurrently; (iii) a wide, HBM-like electrical interface so that the stack can sit on the same interposer as HBM; and (iv) a base die that hosts ECC, address translation and buffering. The public value proposition is HBM-class read bandwidth with 8–16× the capacity per stack — on the order of 512 GB for a first-generation stack — at similar cost, while accepting that write bandwidth and endurance are far lower than DRAM. Because specifications are not yet final, this paper treats every HBF parameter as an explicit, swept assumption (Table 2).");
H2("2.4 Related work");
P("**Weight offloading.** FlexGen [12] and DeepSpeed-Inference [14] offload weights and KV to host DRAM and SSD, maximising throughput for latency-insensitive batch jobs; *LLM in a Flash* [13] exploits activation sparsity to run models larger than DRAM from smartphone flash. These systems traverse PCIe-class links (tens of GB/s) and thus cannot meet interactive TPOT targets for large models. Cambricon-LLM [31] couples an NPU with a flash chiplet with in-flash computation for on-device 70B inference, and in-storage attention offloading [32] moves attention into the SSD. HBF differs by placing flash *on the accelerator package* with two orders of magnitude more bandwidth.");
P("**Serving systems.** PagedAttention/vLLM [6] manages the KV cache like virtual memory; Orca [7] introduced iteration-level scheduling; Sarathi-Serve [8] chunks prefills to bound decode stalls; DistServe [9], Splitwise [10] and Mooncake [11] disaggregate prefill and decode and treat the KV cache as a first-class, transferable object. **KV and weight compression** — multi-query and grouped-query attention [19], [20], multi-head latent attention [17], weight quantisation [39], [40] — reduce bytes moved and are complementary to HBF. **Memory disaggregation** through CXL [33], [34] extends capacity at DRAM latency but at link bandwidth well below HBM. **Processing near or in memory** [35], [36], including in-flash bulk operations [42], reduces data movement by computing at the memory. **SRAM-centric dataflow** machines [43] avoid off-chip memory on the critical path by distributing weights across many chips. To our knowledge no prior work has analysed on-package flash jointly with the prefill/decode split, MoE batch effects, endurance and disaggregated serving; this paper fills that gap.");

// =====================================================================
H1("3. Prefill and Decode Mechanics on a Tiered HBM–HBF Package");
P("We model an accelerator package with effective compute throughput Π (FLOP/s), an HBM tier of capacity *C*_{M} and bandwidth β_{M}, and an HBF tier of capacity *C*_{F} and read bandwidth β_{F}. Model weights occupy *W* bytes (*P* parameters at *b* bytes each); a token's KV footprint is κ bytes. Each tier is accessed through independent channels, so transfers from the two tiers overlap.");
H2("3.1 Prefill");
P("For a prompt chunk of *L* tokens processed against weights resident in tier *X*, the linear-layer FLOPs are approximately 2*P*_{act}*L* and the weight bytes read are *W*_{touched}(*L*). The operational intensity with respect to weight traffic is therefore");
EQ("I_{prefill} = 2 P_{act} L / W_{touched}(L) ≈ 2L / b   (dense models)", 1);
P("which for FP8 weights (*b* = 1) equals 2*L* FLOP/byte. Prefill on tier *X* is compute-bound whenever *I* exceeds the ridge point Π/β_{X}. Solving for *L* gives the **minimum compute-bound chunk size**");
EQ("L^{*} = Π · W_{touched} / (2 β_{X} P_{act})", 2);
P("For a dense model *W*_{touched} = *W* and *L*^{*} = Π*b*/(2β_{X}), independent of model size. For an MoE model, a sufficiently long chunk touches every expert, so *W*_{touched} → *W*_{total} while *P*_{act} remains small; *L*^{*} is therefore inflated by the sparsity ratio *P*/*P*_{act}. This is the first mechanism by which MoE architectures and flash bandwidth interact: an MoE prefill on HBF needs thousands of tokens per weight pass to remain compute-bound (Section 6.4).");
H2("3.2 Decode");
P("One decode iteration over a batch of *B* sequences with mean context *S* reads the touched weights once and each sequence's KV cache once. With weights in HBF and KV in HBM, the step time is");
EQ("t_{step} = max( W_{touched}(B)/β_{F} ,  B·S·κ/β_{M} ,  F(B,S)/Π ) + t_{comm}", 3);
P("where *F*(*B*,*S*) = *B*[2*P*_{act} + 4*S*·*d*_{attn}·*n*_{layers}] counts linear and attention FLOPs and *t*_{comm} captures collective latency when the model is sharded across *N* packages (we use 2*n*_{layers}·α with α = 10 µs). In an HBM-only package both weight and KV traffic share β_{M}, so the first two terms add rather than overlap. Throughput per package is *B*/(*N*·*t*_{step}), and the TPOT seen by each user is *t*_{step}. The batch is bounded by capacity:");
EQ("B·S·κ / N ≤ η C_{M},     W / N ≤ η C_{F}   (hybrid);     (W + B·S·κ)/N ≤ η C_{M}   (HBM-only)", 4);
P("with usable fraction η = 0.9. Equation (4) exposes the central trade-off: moving weights to HBF frees HBM capacity for KV — enlarging the feasible batch — while Equation (3) shows that it also moves weight traffic onto the slower flash channel.");
H2("3.3 Expert-touch probability in MoE decode");
P("In an MoE layer with *E* experts and top-*k* routing, a single token reads only *k* experts. A batch of *B* tokens under approximately uniform routing touches each expert with probability");
EQ("p_{touch}(B) = 1 − (1 − k/E)^{B}", 5);
P("so the expert bytes read per layer scale as *p*_{touch}(*B*)·*W*_{experts}. Fig. 3(b) plots Equation (5): for *E* = 384, *k* = 8 (a configuration representative of current trillion-parameter open models [44]), a batch of 32 already touches 49% of experts and a batch of 128 touches 93%. Consequently the decode weight traffic of a large MoE model is *not* a fixed 2*P*_{act} bytes per step: it grows almost linearly with batch until it saturates at the full expert set. On a flash tier with lower bandwidth than HBM, this traffic growth, rather than KV, becomes the binding constraint (Section 6.2).");
H2("3.4 Latency hiding and required concurrency");
P("Flash read latency is three orders of magnitude above DRAM. By Little's law [38], sustaining bandwidth β with an access latency *t*_{R} requires *Q* outstanding reads of size *s*:");
EQ("Q = β · t_{R} / s", 6);
P("For β = 0.8 TB/s per stack, *t*_{R} = 20 µs and *s* = 16 KiB, *Q* ≈ 977 concurrent page reads per stack (Fig. 6). This is the deepest architectural requirement of HBF: it is only an HBM-class device if the array is partitioned into on the order of a thousand independently operating read units, and if the base-die controller and the accelerator's DMA engines can keep that many requests in flight.");
P("For dense models the weight access stream is entirely deterministic: the weights of layer *ℓ*+1 can be prefetched while layer *ℓ* computes. A decode step of 30–50 ms spread over 60–80 layers allots ≈0.5 ms per layer, which is 10–25× larger than *t*_{R}; a one-layer-ahead prefetch therefore hides flash latency completely. MoE layers break this determinism: which experts to fetch is known only after the router executes. If expert fetches are issued only then, each layer exposes one *t*_{R}, adding *n*_{layers}·*t*_{R} ≈ 61 × 20 µs ≈ 1.2 ms (or 3.1 ms at 50 µs) per token — 2.4–6% of a 50 ms TPOT budget. Speculative expert prefetch based on the previous layer's routing, or fetching all experts once the batch is large enough that *p*_{touch} ≈ 1, removes this penalty.");

// =====================================================================
H1("4. Integrating HBF into Accelerator Systems");
FIG("figures/fig1_architecture.png", "Reference integration architecture. HBM and HBF stacks share the 2.5D interposer with the accelerator die, which integrates separate PHYs and a flash-aware DMA/prefetch engine. HBM holds the mutable KV cache and activations; HBF holds read-mostly model weights and write-once prefix caches. Host DRAM, CXL pools and SSDs form off-package tiers for spill, warm prefixes and checkpoints.");
H2("4.1 Physical integration");
P("HBF is designed to share the interposer and, ideally, the PHY footprint of HBM, so that a package designer can trade HBM sites for HBF sites without redesigning the substrate. Three physical issues require attention. **Beachfront and stack height:** NAND core dies are thinner per bit than DRAM but a high-capacity stack may contain more dies; co-planarity with HBM stacks and the accelerator die constrains die thinning and TSV pitch. **Thermals:** NAND data retention degrades and raw bit-error rate rises with temperature, while the accelerator die adjacent to the stack can dissipate 1 kW. Placing HBF on the cooler edge of the interposer, adding thermal isolation or increasing ECC strength and refresh frequency are the available mitigations; this is an under-studied co-design problem. **Power delivery:** NAND programming requires high internal voltages generated by charge pumps; during bulk weight loading the program current of thousands of concurrently active planes must be budgeted against the package power envelope.");
H2("4.2 Controller and interface");
P("A conventional SSD controller hides NAND behind a block interface with a page-mapped flash translation layer (FTL), garbage collection and deep request queues. HBF's base die must instead deliver deterministic, low-overhead *streaming reads* of large, aligned objects. We argue for a lightweight **read-optimised controller** with the following properties:");
BUL([
  "**Region-mapped translation.** Weights are written once in large contiguous extents; a coarse region table (megabyte granularity) replaces the page-level FTL, eliminating mapping-table DRAM on the base die.",
  "**Load–seal–read life-cycle.** Regions are bulk-programmed when a model is loaded, sealed read-only and erased only when the model is replaced, so garbage collection and write amplification are near zero for the weight tier.",
  "**On-die LDPC with soft-decision fallback.** Hard-decision decoding at line rate on the base die keeps the common-case latency low; read-retry and soft decoding are exceptions handled out of band.",
  "**Read-disturb and retention management.** Weights are read billions of times; the controller must track per-block read counts and transparently relocate blocks that approach disturb limits — a background *write* source that must be included in the endurance budget.",
  "**Descriptor-based bulk DMA.** The accelerator issues descriptors of the form ⟨region, offset, length, destination⟩; the controller expands them into thousands of plane-level reads, satisfying Equation (6).",
]);
H2("4.3 Data placement");
P("Table 1 assigns each inference data object to a tier according to its access pattern, write rate, and latency sensitivity. The governing principle, which Section 6.5 quantifies, is that **HBF should hold only data that is written rarely and read in large, predictable streams**.");
TBL("Placement of inference data objects in a tiered HBM–HBF package.",
  ["Data object", "Access pattern", "Write frequency", "Preferred tier", "Rationale"],
  [
    ["Dense / shared weights", "Sequential stream, every step", "Per model update", "HBF (or HBM if it fits)", "Read-only; deterministic prefetch hides tR"],
    ["MoE expert weights", "Router-dependent, batch-dependent", "Per model update", "HBF; hot experts cached in HBM", "Largest capacity consumer; sparsity tolerates flash bandwidth"],
    ["KV cache (active)", "Read every step, appended every token", "Continuous", "HBM", "Mutable; would exhaust NAND endurance in days"],
    ["Prefix / system-prompt KV", "Read many times", "Written once per prefix", "HBF or CXL", "Write-once-read-many; capacity-hungry"],
    ["LoRA / adapters", "Per-request selection", "Rare", "HBF, staged to HBM", "Many variants, small each"],
    ["Activations, workspace", "Random, short-lived", "Every layer", "SRAM / HBM", "Latency-critical and write-heavy"],
    ["Checkpoints, cold KV", "Bulk", "Occasional", "SSD / object store", "Off-package; capacity only"],
  ], [1700, 1900, 1400, 1800, 2226]);
H2("4.4 Software stack");
P("Integrating HBF changes four layers of the serving software. (1) **Model loader:** checkpoints must be laid out in HBF in the exact order the kernels consume them, with per-layer extents aligned to plane boundaries. (2) **Memory allocator:** the runtime must manage two physical address spaces with different granularity and mutability — a paged KV allocator [6] in HBM and an immutable region allocator in HBF. (3) **Compiler and kernels:** prefetch of layer *ℓ*+1 must be scheduled explicitly and overlapped with compute; for MoE, a speculative expert prefetcher keyed on router history is required. (4) **Scheduler:** the batch-size controller must be aware that in MoE decode the step time grows with *B* through Equation (5), not only through KV traffic, and must therefore choose *B* against the TPOT objective using a cost model of the form of Equation (3).");
H2("4.5 Reliability and endurance budget");
P("The endurance-limited lifetime of an HBF stack of capacity *C* with rated P/E cycles *N*_{PE}, write amplification factor WAF and sustained write rate *r* is");
EQ("T_{life} = N_{PE} · C / (WAF · r)", 7);
P("Section 6.5 evaluates Equation (7) for representative write sources. Because weights are written only when a model is deployed or updated, while the KV cache is written on every token, the two differ in write rate by three to four orders of magnitude — the quantitative basis of the placement in Table 1.");

// =====================================================================
H1("5. Evaluation Methodology");
P("Because HBF silicon is not yet generally available, we evaluate it with an analytical model that implements Equations (2)–(7). Analytical roofline-style models are standard in early-stage architecture studies [4], [5] and, crucially, make every assumption explicit and reproducible. The model and figure scripts are released with this paper (see Data Availability). Table 2 lists device assumptions and Table 3 the workloads.");
TBL("Device and package assumptions (per package of eight memory sites).",
  ["Parameter", "Value", "Notes"],
  [
    ["Effective compute Π", "1.5 PFLOP/s (FP8)", "2.5 PFLOP/s peak × 0.6 utilisation"],
    ["HBM stack", "36 GB, 1.0 TB/s", "12-high HBM3E class"],
    ["HBF stack (baseline)", "512 GB, 0.8 TB/s read", "Projection; 80% of HBM bandwidth"],
    ["HBF stack (optimistic)", "512 GB, 1.6 TB/s read", "Sensitivity study, Section 6.3"],
    ["Configurations", "HBM-8, HBM-6/HBF-2, HBM-4/HBF-4, HBM-2/HBF-6", "Stacks of each type"],
    ["Usable capacity η", "0.9", "Remainder: runtime, fragmentation"],
    ["Collective latency α", "10 µs per collective", "Two collectives per layer when N > 1"],
    ["NAND read latency tR", "3–50 µs", "Swept (Fig. 6)"],
    ["P/E endurance", "3k (TLC), 30k (pSLC)", "pSLC uses one-third capacity; WAF = 2"],
  ], [2600, 3300, 3126]);
TBL("Workload models.",
  ["Model", "Total params", "Active params", "Layers", "Experts (E, k)", "KV bytes/token"],
  [
    ["Dense-70B", "70 B", "70 B", "80", "—", "160 KiB (GQA-8, FP8)"],
    ["MoE-1T", "1.0 T", "32 B", "61", "384, 8", "68.6 KiB (MLA, BF16)"],
  ], [1600, 1400, 1400, 1000, 1500, 2126]);
P("Weights are stored in FP8. For each (model, configuration, *N*, *S*) we sweep the decode batch *B* and report the maximum throughput per package whose step time satisfies TPOT ≤ 50 ms; infeasible combinations (capacity violated at every *B*, or SLO violated at *B* = 1) are reported as such. Tensor/expert parallelism divides weight, KV and FLOP work evenly across *N* packages and adds *t*_{comm}; this idealisation favours scale-out and therefore makes the comparison conservative with respect to HBF. **Limitations:** the model omits energy, interconnect bandwidth contention, routing imbalance (Equation (5) assumes uniform routing; real routers are skewed, which reduces *p*_{touch} and favours HBF), ECC tail latency and queueing; it should be read as a bound-and-trend analysis to be validated against silicon when available.");

// =====================================================================
H1("6. Results");
H2("6.1 Roofline view");
FIG("figures/fig2_roofline.png", "Roofline of a 1.5 PFLOP/s package with an 8 TB/s HBM ceiling and 3.2/4.8 TB/s HBF ceilings (four or six HBF stacks). Decode at small batch lies deep in the bandwidth-bound region; prefill with chunks of a few hundred tokens or more is compute-bound even on the HBF ceiling.", 5.6);
P("Fig. 2 places prefill and decode on a roofline. The ridge point moves from 188 FLOP/byte for an all-HBM package to 469 FLOP/byte when weights stream from four HBF stacks. Decode with *B* ≤ 32 operates at 2–64 FLOP/byte and is bandwidth-bound on any tier, so every byte/s of weight bandwidth translates directly into tokens/s. Dense prefill, at 2*L* FLOP/byte, crosses the HBF ridge for *L* ≥ 235 and is thereafter insensitive to the memory technology. The roofline thus predicts the qualitative split developed below: **HBF is almost free for prefill and expensive for latency-bound decode**.");
H2("6.2 Decode under a latency objective");
FIG("figures/fig3_moe_decode.png", "(a) TPOT versus decode batch for MoE-1T on N = 4 packages at 32k context. The HBM-8 curve terminates at B = 14 because weights plus KV exhaust HBM; hybrid configurations sustain larger batches but their step time rises with B because more experts are touched. (b) Expert-touch probability, Equation (5), for three MoE geometries.");
P("**Dense-70B.** A 70B model in FP8 fits in a single HBM-8 package, and HBF never helps. At 8k context and *N* = 1, HBM-8 delivers 4,308 tokens/s per package within the SLO versus 2,980 for HBM-4/HBF-4 and 1,490 for HBM-2/HBF-6; the pattern persists for every *N* and at 32k context (1,077 vs. 745 vs. 373). Replacing HBM stacks with HBF removes KV capacity and KV bandwidth that a model of this size needs more than it needs weight capacity. This negative result is important for the debate: **HBF is not a general-purpose HBM substitute**.");
P("**MoE-1T.** The picture changes for a trillion-parameter MoE (Fig. 3, Fig. 4, Table 4). An HBM-only deployment is infeasible below four packages. With HBF the model fits in a *single* package, which serves 152 (HBM-4/HBF-4) to 248 (HBM-2/HBF-6) tokens/s at 46–49 ms TPOT. At four packages and 32k context, HBM-8 is capacity-starved — after storing 250 GB of weights per package only ≈9 GB remains for KV, limiting the batch to 12 — and delivers 321 tokens/s per package, while HBM-2/HBF-6 sustains a batch of 108 and 562 tokens/s per package, a **1.75× advantage**. At 8k context the two are at parity (574 vs. 562) because shorter contexts relieve HBM capacity. From eight packages upward, the aggregate HBM capacity of the HBM-only design is sufficient and its higher bandwidth wins (1,645 vs. 1,307 tokens/s per package at 32k).");
FIG("figures/fig4_slo_throughput.png", "Best decode throughput per package subject to TPOT ≤ 50 ms for MoE-1T versus the number of packages N, at 8k and 32k context. Missing bars denote infeasible configurations. HBF configurations are the only feasible options at N ≤ 2 and lead at N = 4 for long contexts; HBM-only leads once N is large enough to hold weights plus KV in HBM.");
TBL("Decode throughput per package (tokens/s) for MoE-1T at TPOT ≤ 50 ms, 32k context. Batch size in parentheses; the binding resource is the tier listed in the last column for the best hybrid entry.",
  ["Packages N", "HBM-8", "HBM-6/HBF-2", "HBM-4/HBF-4", "HBM-2/HBF-6", "Binding resource"],
  [
    ["1", "infeasible", "infeasible", "152 (7)", "248 (12)", "HBF weight bandwidth"],
    ["2", "infeasible", "74 (7)", "169 (16)", "287 (27)", "HBF weight bandwidth"],
    ["4", "321 (12)", "85 (16)", "223 (43)", "562 (108)", "HBF weight bandwidth"],
    ["8", "1,645 (421)", "112 (43)", "1,307 (421)", "835 (214)", "HBF / HBM-KV"],
    ["16", "2,500 (1,289)", "2,000 (1,289)", "1,668 (824)", "835 (421)", "HBM-KV bandwidth"],
  ], [1200, 1500, 1500, 1500, 1500, 1826]);
P("Two further observations matter for designers. First, **HBF stack count should be sized by bandwidth, not by capacity**: HBM-6/HBF-2 provides 1 TB of flash — enough to hold the model — but its 1.6 TB/s of weight bandwidth makes it the worst configuration at every *N* ≤ 8. Second, the binding resource in every capacity-constrained hybrid case is HBF weight bandwidth driven by the expert-touch growth of Equation (5), not KV traffic. Techniques that reduce *p*_{touch} — expert-affinity batching, hot-expert caching in HBM, or routing-aware request grouping — therefore have outsized value on HBF systems.");
H2("6.3 Sensitivity to HBF bandwidth");
P("Doubling HBF read bandwidth to 1.6 TB/s per stack — plausible for later generations with more sub-array parallelism or a wider interface — raises single-package MoE-1T throughput from 152 to 348 tokens/s (HBM-4/HBF-4) and from 248 to 589 tokens/s (HBM-2/HBF-6). At *N* = 4 and 32k context, HBM-4/HBF-4 reaches 1,342 tokens/s per package with a batch of 214, **4.2×** the HBM-8 result, and HBM-2/HBF-6 becomes HBM-KV-bound at 836 tokens/s. The benefit of HBF is thus highly elastic in its bandwidth: at parity with HBM bandwidth the hybrid would dominate every capacity-constrained regime, whereas at half of HBM bandwidth its advantage would vanish except for feasibility at small *N*.");
H2("6.4 Prefill");
TBL("Minimum compute-bound prefill chunk L* from Equation (2) and time to prefill an 8,192-token prompt (compute-limited and weight-traffic-limited components).",
  ["Model / configuration", "L* (tokens)", "t_compute (s)", "t_weights (s)", "Bound at L = 8,192"],
  [
    ["Dense-70B / HBM-8", "94", "0.82", "0.009", "Compute"],
    ["Dense-70B / HBM-4/HBF-4", "234", "0.82", "0.022", "Compute"],
    ["MoE-1T / HBM-8 (N ≥ 4)", "2,930", "0.44", "0.125", "Compute"],
    ["MoE-1T / HBM-4/HBF-4", "7,324", "0.44", "0.313", "Compute (marginal)"],
    ["MoE-1T / HBM-2/HBF-6", "4,883", "0.44", "0.208", "Compute"],
  ], [2800, 1400, 1500, 1500, 1826]);
P("Table 5 confirms that dense prefill is compute-bound on HBF for any chunk above a few hundred tokens. MoE prefill is more demanding: with four HBF stacks a chunk must contain at least 7,324 tokens for flash weight traffic to be hidden behind computation, because a long chunk touches every expert while only 32B parameters are active per token. Chunked-prefill schedulers [8] that use 512–2,048-token chunks to bound decode interference would, on an HBF-weight package serving a large MoE, become weight-bandwidth-bound. This creates a direct tension between chunked prefill and HBF that **prefill–decode disaggregation resolves**: a dedicated prefill node can batch prompts into large chunks without harming any decode TPOT.");
H2("6.5 Endurance");
FIG("figures/fig5_endurance.png", "Endurance-limited lifetime of a 512 GB HBF stack versus sustained write rate (WAF = 2). Vertical markers indicate a daily full rewrite of a 256 GB weight shard, sustained KV writes (≈1 GB/s) and KV swap-in bursts (≈20 GB/s).", 5.6);
P("Fig. 5 evaluates Equation (7). Rewriting a stack's 256 GB weight shard once per day (0.003 GB/s) consumes the endurance of a TLC stack in 2,963 days (8.1 years) and of a pSLC configuration in 27 years — comfortably beyond a five-year service life even with daily model updates. In contrast, a sustained KV write rate of 1 GB/s — typical of a package generating several thousand tokens/s with 70–160 KiB of KV per token — wears out TLC in **8.9 days** and pSLC in 30 days; KV swap-in bursts at 20 GB/s would do so in hours. The KV cache must therefore never be placed in HBF as a general-purpose spill target. A write-once prefix cache — a system prompt or shared document reused by many requests [11] — has a write rate set by the arrival of *new* prefixes and is acceptable when that rate is below ≈20 MB/s per stack.");
H2("6.6 Concurrency requirements");
FIG("figures/fig6_littles_law.png", "Concurrent reads in flight per stack needed to sustain 0.8 TB/s, from Equation (6), for 4, 16 and 64 KiB read units.", 5.6);
P("Fig. 6 translates the bandwidth target into a parallelism requirement. With 16 KiB pages, a stack must keep 146 reads in flight at *t*_{R} = 3 µs (low-latency SLC-class arrays), 977 at 20 µs and 2,441 at 50 µs; with 4 KiB units these numbers quadruple to 586–9,766. A 16-die stack with four planes per die offers only 64 conventional planes, more than an order of magnitude short of the 20 µs requirement. Hence HBF requires either many more independently sensed sub-arrays per die, larger effective read units through multi-plane reads, lower-latency array designs, or all three — a central VLSI and circuit-design challenge for the technology.");

// =====================================================================
H1("7. The System Architecture Debate");
P("The inference memory hierarchy is the subject of an active debate between at least five architectural positions. Table 6 summarises them; we then argue where each is strongest, using the results of Section 6.");
TBL("Competing architectural positions for the LLM inference memory hierarchy.",
  ["Position", "Core idea", "Strengths", "Weaknesses", "Best regime"],
  [
    ["A. HBM scale-out", "More HBM per package (HBM4/4E), shard across many packages", "Highest bandwidth; mature software", "Capacity cost; packages sized by memory, not compute", "Latency-critical serving of models that fit"],
    ["B. On-package HBF", "Flash stacks beside HBM on interposer", "8–16× capacity/stack; fewer packages", "Read-only practicality; tR; endurance; immature", "Large MoE, small N, prefill nodes"],
    ["C. CXL / host pools", "Off-package DRAM or flash over CXL/PCIe", "Elastic, shareable capacity", "Link bandwidth ≪ HBM", "KV spill, prefix caches, batch jobs"],
    ["D. Near/in-memory compute", "Compute inside HBM or flash base dies", "Avoids moving weights; high internal BW", "Programmability; thermal; vendor lock-in", "Attention, GEMV, expert FFN offload"],
    ["E. SRAM-centric dataflow", "Distribute weights across many chips' SRAM", "Lowest latency per token", "Huge chip count for TB-scale models", "Low-latency serving of mid-size models"],
  ], [1500, 1900, 1800, 1900, 1926]);
P("**Position A versus B.** Our results support neither extreme. For models that fit in HBM, Position A is strictly better (Section 6.2, Dense-70B). For trillion-parameter MoE models at small-to-moderate package counts — the regime of on-premises, sovereign, enterprise and edge-datacentre deployments, where package count is fixed by budget or power — Position B is the only option at *N* ≤ 2 and leads by 1.75–4.2× at *N* = 4 for long contexts. At hyperscale, where operators can deploy as many packages as needed, Position A remains competitive on throughput per package, and the argument for HBF shifts to cost per token, which depends on HBF pricing not yet public.");
P("**Position C.** CXL-attached memory [33], [34] offers capacity elasticity and sharing across packages but at link bandwidths one to two orders of magnitude below HBM. Our decode analysis implies it cannot hold per-step weight traffic for interactive serving; it is, however, an excellent home for warm prefix caches and preempted KV, and is complementary to HBF.");
P("**Position D.** Processing-near-memory [35], [36] and in-flash computation [31], [42] would invert the trade-off identified in Section 6.2: if the expert FFN computation for a token were performed in the HBF base die, only activations — kilobytes rather than gigabytes — would cross the interface, removing HBF bandwidth as the binding constraint. The base die that HBF already requires is a natural host for such logic, but thermal limits next to NAND and the need for a programmable, standardised interface are open problems.");
P("**Position E.** SRAM-centric designs [43] minimise TPOT but require a number of chips proportional to model size; for trillion-parameter models the chip count, and thus capital cost, is prohibitive. They and HBF sit at opposite corners of the capacity–latency plane.");
H2("7.1 A heterogeneous fleet: HBF for prefill, HBM for decode");
P("Combining Sections 6.2 and 6.4 yields a concrete architectural proposal. Prefill is compute-bound when chunks are large, needs little *persistent* KV capacity (KV is shipped to decode nodes immediately [9]–[11]), and tolerates HBF weight bandwidth once *L* ≥ *L*^{*}. Decode is bandwidth- and KV-capacity-bound and is hurt by moving weight traffic to flash unless capacity, not bandwidth, is the binding constraint. We therefore propose that **prefill pools be built from HBF-heavy packages** (e.g., HBM-2/HBF-6), which can hold an entire trillion-parameter model per package and thus prefill without any model-parallel communication, while **decode pools use HBM-rich packages**, or HBM-4/HBF-4 hybrids when the deployment size is fixed and small. This assignment uses each memory technology where its weaknesses are masked by the phase's characteristics, and it is compatible with KV-centric serving architectures already in production [11].");

// =====================================================================
H1("8. Open Research Challenges");
BUL([
  "**Array-level parallelism (VLSI).** Achieving ≈10^{3} concurrent reads per stack (Section 6.6) at acceptable area, power and peripheral-circuit overhead; low-latency sensing schemes and sub-array partitioning in 3D NAND.",
  "**Thermal co-design.** Retention, bit-error rate and ECC strength of NAND operated next to a kilowatt-class accelerator; interposer floor-planning and thermal isolation.",
  "**Interface standardisation.** An HBM-compatible PHY with a read-optimised command set (bulk descriptors, sealed regions) to allow mixed HBM/HBF packages from multiple vendors [30].",
  "**Expert-aware scheduling.** Batch formation that minimises *p*_{touch} (Equation (5)) through routing-affinity grouping; hot-expert caching policies between HBM and HBF.",
  "**Speculative expert prefetch.** Predicting next-layer expert selection to hide *t*_{R} in MoE decode without wasting flash bandwidth.",
  "**Reliability under read-heavy stress.** Read-disturb characterisation at billions of reads per block, and the background relocation write budget it implies.",
  "**Compute in the base die.** Programmable near-flash GEMV/FFN engines that remove weight traffic from the interface (Position D).",
  "**Validation.** Cycle-accurate and silicon-level validation of the analytical model presented here, including energy per token, which this study does not quantify.",
]);

// =====================================================================
H1("9. Conclusion");
P("High Bandwidth Flash reframes the LLM inference memory problem from \"how much HBM can be afforded\" to \"which data deserves DRAM\". By deriving and evaluating the governing equations for prefill and decode on a tiered HBM–HBF package, this paper has shown that HBF is not a drop-in substitute for HBM: it lowers per-package decode throughput for models that already fit in HBM, and its value collapses if its read bandwidth falls well below HBM's or if it is used for mutable data. Where it is used for what flash does well — holding enormous, read-mostly weight sets — it changes what is feasible: a trillion-parameter MoE fits in a single package, capacity-constrained deployments gain 1.75–4.2× throughput per package, prefill nodes can be built without model parallelism, and daily model updates consume only a small fraction of endurance. The decisive technical questions are now at the device and circuit level — sub-array parallelism, latency and thermal robustness — and at the system level — expert-aware scheduling and prefill–decode disaggregation. We hope the analytical framework and open model presented here help both communities reason about them quantitatively.");

// ---------- declarations ----------
H1("Declarations");
P("**Funding.** This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.", { noIndent: true });
P("**Conflict of interest.** The author declares no competing interests.", { noIndent: true });
P("**Data and code availability.** The analytical model (hbf_model.py), figure-generation script (make_figures.py) and raw results (results.json) are available from the author's repository alongside this manuscript.", { noIndent: true });
P("**Author contribution.** M. Nisha Angeline: conceptualisation, methodology, analysis, visualisation, writing — original draft, review and editing.", { noIndent: true });

// ---------- references ----------
H1("References");
const refs = [
  "A. Vaswani *et al.*, \"Attention is all you need,\" in *Proc. Adv. Neural Inf. Process. Syst. (NeurIPS)*, 2017, pp. 5998–6008.",
  "J. Kaplan *et al.*, \"Scaling laws for neural language models,\" arXiv:2001.08361, 2020.",
  "A. Gholami, Z. Yao, S. Kim, C. Hooper, M. W. Mahoney, and K. Keutzer, \"AI and memory wall,\" *IEEE Micro*, vol. 44, no. 3, pp. 33–39, 2024.",
  "R. Pope *et al.*, \"Efficiently scaling transformer inference,\" in *Proc. Mach. Learn. Syst. (MLSys)*, 2023.",
  "S. Williams, A. Waterman, and D. Patterson, \"Roofline: An insightful visual performance model for multicore architectures,\" *Commun. ACM*, vol. 52, no. 4, pp. 65–76, 2009.",
  "W. Kwon *et al.*, \"Efficient memory management for large language model serving with PagedAttention,\" in *Proc. ACM Symp. Operating Syst. Principles (SOSP)*, 2023.",
  "G.-I. Yu, J. S. Jeong, G.-W. Kim, S. Kim, and B.-G. Chun, \"Orca: A distributed serving system for transformer-based generative models,\" in *Proc. USENIX OSDI*, 2022.",
  "A. Agrawal *et al.*, \"Taming throughput–latency tradeoff in LLM inference with Sarathi-Serve,\" in *Proc. USENIX OSDI*, 2024.",
  "Y. Zhong *et al.*, \"DistServe: Disaggregating prefill and decoding for goodput-optimized large language model serving,\" in *Proc. USENIX OSDI*, 2024.",
  "P. Patel *et al.*, \"Splitwise: Efficient generative LLM inference using phase splitting,\" in *Proc. ACM/IEEE Int. Symp. Comput. Archit. (ISCA)*, 2024.",
  "R. Qin *et al.*, \"Mooncake: A KVCache-centric disaggregated architecture for LLM serving,\" in *Proc. USENIX Conf. File Storage Technol. (FAST)*, 2025.",
  "Y. Sheng *et al.*, \"FlexGen: High-throughput generative inference of large language models with a single GPU,\" in *Proc. Int. Conf. Mach. Learn. (ICML)*, 2023.",
  "K. Alizadeh *et al.*, \"LLM in a flash: Efficient large language model inference with limited memory,\" in *Proc. Annu. Meeting Assoc. Comput. Linguistics (ACL)*, 2024.",
  "R. Y. Aminabadi *et al.*, \"DeepSpeed-Inference: Enabling efficient inference of transformer models at unprecedented scale,\" in *Proc. Int. Conf. High Perform. Comput., Netw., Storage Anal. (SC)*, 2022.",
  "N. Shazeer *et al.*, \"Outrageously large neural networks: The sparsely-gated mixture-of-experts layer,\" in *Proc. Int. Conf. Learn. Represent. (ICLR)*, 2017.",
  "W. Fedus, B. Zoph, and N. Shazeer, \"Switch Transformers: Scaling to trillion parameter models with simple and efficient sparsity,\" *J. Mach. Learn. Res.*, vol. 23, no. 120, pp. 1–39, 2022.",
  "DeepSeek-AI, \"DeepSeek-V2: A strong, economical, and efficient mixture-of-experts language model,\" arXiv:2405.04434, 2024.",
  "DeepSeek-AI, \"DeepSeek-V3 technical report,\" arXiv:2412.19437, 2024.",
  "J. Ainslie *et al.*, \"GQA: Training generalized multi-query transformer models from multi-head checkpoints,\" in *Proc. Conf. Empirical Methods Natural Lang. Process. (EMNLP)*, 2023.",
  "N. Shazeer, \"Fast transformer decoding: One write-head is all you need,\" arXiv:1911.02150, 2019.",
  "T. Dao, D. Y. Fu, S. Ermon, A. Rudra, and C. Ré, \"FlashAttention: Fast and memory-efficient exact attention with IO-awareness,\" in *Proc. NeurIPS*, 2022.",
  "Y. Leviathan, M. Kalman, and Y. Matias, \"Fast inference from transformers via speculative decoding,\" in *Proc. ICML*, 2023.",
  "JEDEC Solid State Technology Association, *High Bandwidth Memory DRAM (HBM3)*, JESD238, 2022.",
  "D. U. Lee *et al.*, \"A 1.2 V 8 Gb 8-channel 128 GB/s high-bandwidth memory (HBM) stacked DRAM with effective microbump I/O test methods using 29 nm process and TSV,\" in *IEEE Int. Solid-State Circuits Conf. (ISSCC) Dig. Tech. Papers*, 2014, pp. 432–433.",
  "M. O'Connor *et al.*, \"Fine-grained DRAM: Energy-efficient DRAM for extreme bandwidth systems,\" in *Proc. IEEE/ACM Int. Symp. Microarchitecture (MICRO)*, 2017.",
  "R. Micheloni, L. Crippa, and A. Marelli, *Inside NAND Flash Memories*. Dordrecht, The Netherlands: Springer, 2010.",
  "Y. Cai, S. Ghose, E. F. Haratsch, Y. Luo, and O. Mutlu, \"Error characterization, mitigation, and recovery in flash-memory-based solid-state drives,\" *Proc. IEEE*, vol. 105, no. 9, pp. 1666–1704, 2017.",
  "H. Tanaka *et al.*, \"Bit cost scalable technology with punch and plug process for ultra high density flash memory,\" in *Proc. IEEE Symp. VLSI Technol.*, 2007, pp. 14–15.",
  "Sandisk Corporation, \"Sandisk unveils High Bandwidth Flash (HBF) memory technology for AI inference,\" Investor Day press release, Feb. 2025.",
  "Sandisk Corporation and SK hynix, \"Memorandum of understanding to standardize High Bandwidth Flash (HBF),\" press release, Aug. 2025.",
  "Z. Yu *et al.*, \"Cambricon-LLM: A chiplet-based hybrid architecture for on-device inference of 70B LLM,\" in *Proc. IEEE/ACM MICRO*, 2024.",
  "X. Pan *et al.*, \"InstInfer: In-storage attention offloading for cost-effective long-context LLM inference,\" arXiv preprint, 2024.",
  "H. Li *et al.*, \"Pond: CXL-based memory pooling systems for cloud platforms,\" in *Proc. ACM Int. Conf. Archit. Support Program. Lang. Oper. Syst. (ASPLOS)*, 2023.",
  "D. Das Sharma, R. Blankenship, and D. S. Berger, \"An introduction to the Compute Express Link (CXL) interconnect,\" *ACM Comput. Surv.*, 2024.",
  "O. Mutlu, S. Ghose, J. Gómez-Luna, and R. Ausavarungnirun, \"Processing data where it makes sense: Enabling in-memory computation,\" *Microprocess. Microsyst.*, vol. 67, pp. 28–41, 2019.",
  "S. Lee *et al.*, \"Hardware architecture and software stack for PIM based on commercial DRAM technology,\" in *Proc. ACM/IEEE ISCA*, 2021.",
  "N. P. Jouppi *et al.*, \"TPU v4: An optically reconfigurable supercomputer for machine learning with hardware support for embeddings,\" in *Proc. ACM/IEEE ISCA*, 2023.",
  "J. D. C. Little, \"A proof for the queuing formula: L = λW,\" *Oper. Res.*, vol. 9, no. 3, pp. 383–387, 1961.",
  "E. Frantar, S. Ashkboos, T. Hoefler, and D. Alistarh, \"GPTQ: Accurate post-training quantization for generative pre-trained transformers,\" in *Proc. ICLR*, 2023.",
  "J. Lin *et al.*, \"AWQ: Activation-aware weight quantization for on-device LLM compression and acceleration,\" in *Proc. MLSys*, 2024.",
  "J. Hoffmann *et al.*, \"Training compute-optimal large language models,\" in *Proc. NeurIPS*, 2022.",
  "J. Park *et al.*, \"Flash-Cosmos: In-flash bulk bitwise operations using inherent computation capability of NAND flash memory,\" in *Proc. IEEE/ACM MICRO*, 2022.",
  "D. Abts *et al.*, \"Think fast: A tensor streaming processor (TSP) for accelerating deep learning workloads,\" in *Proc. ACM/IEEE ISCA*, 2020.",
  "Kimi Team, \"Kimi K2: Open agentic intelligence,\" arXiv preprint, 2025.",
];
refs.forEach((r, i) => body.push(new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { after: 60 },
  indent: { left: 480, hanging: 480 }, children: [new TextRun({ text: `[${i + 1}]\t`, size: 19 }), ...runs(r, { size: 19 })],
  tabStops: [{ type: TabStopType.LEFT, position: 480 }] })));

// =====================================================================
const doc = new Document({
  creator: "M. Nisha Angeline",
  title: "High Bandwidth Flash for Large Language Model Inference",
  description: "Prefill–decode mechanics, system integration and architecture debate for HBF",
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT, color: "1F3A5F" }, paragraph: { spacing: { before: 300, after: 140 }, keepNext: true, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 23, bold: true, italics: true, font: FONT }, paragraph: { spacing: { before: 200, after: 100 }, keepNext: true, outlineLevel: 1 } },
    ],
  },
  numbering: { config: [
    { reference: "bul", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
      style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] },
    { reference: "contrib", levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "(%1)", alignment: AlignmentType.LEFT,
      style: { paragraph: { indent: { left: 720, hanging: 400 } } } }] },
  ] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    headers: { default: new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT,
      children: [new TextRun({ text: "High Bandwidth Flash for LLM Inference", italics: true, size: 18, color: "666666" })] })] }) },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: [PageNumber.CURRENT], size: 18 })] })] }) },
    children: body,
  }],
});
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync("HBF_LLM_Inference_Paper.docx", buf); console.log("written"); });
