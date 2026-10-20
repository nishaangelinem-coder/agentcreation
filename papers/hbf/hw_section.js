// Section 7 of the manuscript: near-flash multiplier engines (data-driven from hw/results).
const fs = require("fs");
const T = JSON.parse(fs.readFileSync(__dirname + "/hw/results/hw_tables.json", "utf8"));
for (const k of Object.keys(T)) if (Array.isArray(T[k])) T[k].forEach((r) => { for (const f of ["arch", "fam", "name"]) if (typeof r[f] === "string") r[f] = r[f].replace("a*b", "a×b"); });
const TESTS = fs.readFileSync(__dirname + "/hw/results/functional_tests.txt", "utf8").trim().split("\n");

const f = (x, d = 2) => (x >= 1000 ? Math.round(x).toLocaleString("en-US") : Number(x).toFixed(d));
const g3 = (x) => (x >= 100 ? f(x, 0) : x >= 10 ? f(x, 1) : x >= 1 ? f(x, 2) : f(x, 3));
const row = (tab, d, mode = "area") => T[tab].find((r) => r.design === d && r.mode === mode);
const mac = (d, lib) => T.macs.find((r) => r.design === d && r.lib === lib);
const cap = (d, lib) => T.nmc_capacity.find((r) => r.design === d && r.lib === lib);
const wid = (fam, n, lib) => T.width.find((r) => r.fam === fam && r.n === n && r.lib === lib);
const sysE = (model, B, prefix) => T.sys_energy.find((r) => r.model === model && r.B === B && r.scn.startsWith(prefix));

module.exports = function addHwSection({ P, H1, H2, EQ, FIG, TBL, BUL, NUM }) {
  const s8 = (d) => row("mult8_sky130", d), a8 = (d) => row("mult8_asap7", d);
  const A = a8("bw_array8"), D = a8("bw_dadda8"), DK = a8("bw_dadda8_ks"), Bo = a8("booth8"), BK = a8("booth8_ks"),
    Wl = a8("bw_wallace8"), Be = a8("behav_mult8");
  const As = s8("bw_array8"), Bs = s8("booth8"), BKs = s8("booth8_ks"), DKs = s8("bw_dadda8_ks");
  const best = (key, lib, dir = -1) => {
    const rows = T[`mult8_${lib}`].filter((r) => r.mode === "area" && !r.design.includes("_p"));
    return rows.reduce((b, r) => ((dir < 0 ? r[key] < b[key] : r[key] > b[key]) ? r : b));
  };
  const st = T.struct;
  const dd = (() => { const t = [], e = [];
    for (const lib of ["sky130", "asap7"]) for (const r of T[`mult8_${lib}`].filter((x) => x.mode === "area" && !x.design.includes("_p"))) {
      const q = row(`mult8_${lib}`, r.design, "delay"); t.push((1 - q.tclk / r.tclk) * 100); e.push((q.e_pj / r.e_pj - 1) * 100); }
    return { tmin: Math.min(...t), tmax: Math.max(...t), emax: Math.max(...e) }; })();
  const worst8 = Math.max(...["sky130", "asap7"].flatMap((lib) => ["area", "tclk", "e_pj"].map((k) => row(`mult8_${lib}`, "bw_array8")[k] / best(k, lib)[k])));
  const npass = TESTS.filter((l) => l.includes("PASS")).length;

  H1("7. Near-Flash Compute: Integrating Array Multiplier Engines into the HBF Base Die");
  P("Section 6.2 identified HBF weight bandwidth as the binding resource whenever HBF is used for capacity, and the architecture debate of Section 8 (Position D) argues that computing *inside* the stack would remove weight traffic from the interposer. This section makes that argument concrete at the circuit level. We ask four questions: (i) which multiplier architecture is best suited to a near-flash engine; (ii) how the candidates compare as stand-alone arithmetic units, i.e. *without* HBF; (iii) how many multiply–accumulate (MAC) units an HBF base die can host within realistic area and thermal budgets; and (iv) what the system gains in energy and decode throughput relative to a conventional HBM system and to HBF with far compute.");

  H2("7.1 Placement options");
  TBL("Candidate locations for multiplier hardware in an HBF stack.",
    ["Location", "Mechanism", "Precision", "Main obstacle", "Assessment"],
    [
      ["Logic base die", "Weights stream from page buffers through TSVs into a digital MAC array", "INT4/INT8/FP8, exact", "Area and heat beneath NAND dies", "Feasible; studied here"],
      ["Page-buffer level", "Bit-serial shift/add/bitwise logic in each plane's page buffer [42]", "Low; bit-serial", "Hundreds of cycles per multiply", "Niche (binary/ternary)"],
      ["In NAND array (analog)", "Conductance-weighted string-current summation", "≈4–6 bit, noisy", "ADCs, drift, read disturb", "Not yet viable for LLMs"],
    ], [1500, 2600, 1300, 1800, 1826]);

  H2("7.2 Multiplier architectures and analytical complexity");
  P("All candidates share the same partial-product generation for signed operands, the modified Baugh–Wooley (BW) formulation [45], [60], which avoids sign extension by complementing the partial products that involve exactly one sign bit and adding two constants:");
  EQ("A·B = a_{n−1}b_{n−1}2^{2n−2} + Σ_{i,j<n−1} a_{i}b_{j}2^{i+j} + Σ_{i<n−1} ¬(a_{i}b_{n−1})2^{i+n−1} + Σ_{j<n−1} ¬(a_{n−1}b_{j})2^{j+n−1} + 2^{n} + 2^{2n−1}", 8);
  P("The architectures differ in how they reduce the *n*^{2}-bit partial-product matrix to two rows and in the final carry-propagate adder (CPA). The **array multiplier** accumulates one row per stage with a carry-save adder line, giving a regular, abutting layout but a reduction depth that grows linearly with *n*. **Wallace** [46] and **Dadda** [47] trees reduce all columns concurrently with 3:2 counters; Dadda uses the minimum number of counters per level and the height sequence *d*_{j+1} = ⌊1.5 *d*_{j}⌋ (2, 3, 4, 6, 9, 13, …). **Radix-4 Booth** recoding [48], [49] of *B* halves the number of partial-product rows to *n*/2 at the cost of a recoder and a 5:1 selector per bit; we reduce Booth rows with a Dadda tree and handle sign extension with the same complement-and-constant technique as Equation (8). Each reduction can terminate in a ripple-carry adder (RCA) or a Kogge–Stone (KS) parallel-prefix adder [50]. Table 7 lists the counter counts and logic depths measured on the generated netlists (2-input-gate levels; a full adder counts two levels), which follow the closed forms");
  EQ("Array:  N_{FA} = (n−1)^{2},   N_{HA} = n,   D_{tree} = 2(n−2) + 1", "9a");
  EQ("Dadda:  D_{tree} = 2L + 1,   L = smallest j with d_{j} ≥ n", "9b");
  P("— for example, a 16-bit array needs a 29-level reduction while a 16-bit Dadda tree needs 13 levels, and Booth recoding reduces the 16-bit counter count from 225 to 120 full adders.");
  const sr = (k, lab) => [lab, st[k].fa, st[k].ha, st[k].tree_depth, st[k].depth, st[k].regs || "—"];
  TBL("Analytical structure of the generated multipliers (counts from the netlist generator; depth in 2-input-gate levels).",
    ["Design", "Full adders", "Half adders", "Reduction depth", "Total depth", "Pipeline flops"],
    [
      sr("array4", "Array 4×4"), sr("dadda4", "Dadda 4×4"), sr("booth4", "Booth-4 4×4"),
      sr("array8", "Array 8×8"), sr("wallace8", "Wallace 8×8"), sr("dadda8", "Dadda 8×8"), sr("dadda+ks8", "Dadda+KS 8×8"),
      sr("booth8", "Booth-4 8×8"), sr("booth+ks8", "Booth-4+KS 8×8"),
      sr("array16", "Array 16×16"), sr("wallace16", "Wallace 16×16"), sr("dadda16", "Dadda 16×16"), sr("dadda+ks16", "Dadda+KS 16×16"),
      sr("booth16", "Booth-4 16×16"), sr("booth+ks16", "Booth-4+KS 16×16"),
      sr("array8_p3", "Array 8×8, 3-stage"), sr("dadda8_p2", "Dadda 8×8, 2-stage"), sr("booth8_p2", "Booth 8×8, 2-stage"), sr("array16_p4", "Array 16×16, 4-stage"),
    ], [2300, 1300, 1300, 1400, 1300, 1426]);
  P("With an RCA, the 2*n*-bit carry chain dominates the total depth of every tree multiplier (29 levels at *n* = 8 for array, Wallace and Dadda alike); the tree advantage only becomes visible when it is paired with a parallel-prefix CPA (19 levels for Dadda+KS at *n* = 8). Pipelining cuts the per-stage depth to 15 levels at 8 bits and 31 at 16 bits.");

  H2("7.3 RTL implementation and verification");
  FIG("figures/fig7_nmc_engine.png", "Proposed near-flash MAC engine in the HBF logic base die. Each group of NAND channels delivers decoded weight pages to a systolic MAC tile; activations are broadcast from a small SRAM, and only activations and partial sums cross the interposer to the accelerator.");
  P("All multipliers are produced by an open-source Python generator that emits structural Verilog from the partial-product dot matrix, so that every architecture shares identical partial-product logic and differs only in the reduction and CPA. The generator supports arbitrary widths, signed (BW or Booth) and unsigned operands, and the insertion of pipeline registers after any reduction stage or at the midpoint of the CPA. On top of the multipliers we implement three engines for the base die:");
  BUL([
    "**INT8 MAC processing element (PE):** BW array multiplier, one product register and a 32-bit accumulator, with pass-through registers for the systolic dataflow.",
    "**FP8 (E4M3) [57] MAC with exact accumulation [58]:** a 4×4 unsigned array multiplier on the significands (with the hidden bit), a 5-bit exponent adder, a barrel shifter that aligns the 8-bit significand product to a fixed-point grid with resolution 2^{−18}, and a 48-bit two's-complement accumulator. Because every E4M3 product is an integer multiple of 2^{−18} with magnitude below 2^{18}, Equation (10) is exact and a dot product of up to 2^{11} full-scale terms cannot overflow:",
    "**4×8 output-stationary systolic array:** each column is fed by one NAND channel lane (one weight byte per cycle), each row holds one sequence of the decode batch; PE(*r*, *c*) accumulates Σ_{k} *A*_{rk}*W*_{kc}.",
  ]);
  EQ("w·x = (−1)^{s} · m_{w}m_{x} · 2^{Ew+Ex−20} = (−1)^{s} · [ m_{w}m_{x} ≪ (E_{w}+E_{x}−2) ] · 2^{−18}", 10);
  P("with *s* = *s*_{w} ⊕ *s*_{x}, significands *m* = {*e* ≠ 0, mantissa} and effective exponents *E* = max(*e*, 1).", { noIndent: true });
  P(`Functional verification used Icarus Verilog [56]. Every 4- and 8-bit multiplier was checked **exhaustively** (all 256 or 65,536 operand pairs), 16-bit designs with 300,000 random vectors plus all sign/extreme corner combinations, pipelined variants cycle-accurately against a latency-aligned reference, the FP8 MAC exhaustively for all 65,536 products and on 200 random dot products of length 1–1,024 against an exact rational (Python \`Fraction\`) reference, and the systolic array on 50 random 4×64 by 64×8 matrix products. All ${npass} test suites pass with zero mismatches.`);

  H2("7.4 Synthesis and power methodology");
  P(`Each design was synthesised with Yosys [52] and ABC [53] to two open standard-cell libraries representing the range of plausible base-die processes: **SkyWater 130 nm** [55] (sky130_fd_sc_hd, TT, 1.8 V, 25 °C), standing in for a mature, low-cost logic node, and the **ASAP7** 7 nm predictive library [54] (RVT, TT, 0.7 V), standing in for an advanced node. Every design was mapped twice — *area-driven* (ABC &nf without a delay target) and *delay-driven* (minimum delay target followed by buffering and gate sizing); adders inferred from behavioural code are mapped by Yosys to Brent–Kung parallel-prefix structures [51]. The critical path is ABC's static timing estimate; for registered metrics we add the library flip-flop clock-to-Q and setup times (${f(T.ff_overhead.sky130, 0)} ps for sky130 and ${f(T.ff_overhead.asap7, 0)} ps for ASAP7). Power was obtained from **unit-delay gate-level simulation** of the mapped netlists — which, unlike zero-delay simulation, propagates the spurious transitions (glitches) characteristic of array multipliers — over 3,000 cycles (300 for the 17,000-cell systolic array) for two operand statistics: uniformly random operands and a zero-mean Gaussian distribution with σ = 1/8 of full scale that resembles quantised LLM weights. Switching energy is ½*C*·*V*^{2} per toggle with *C* the sum of fan-out pin capacitances from the liberty file, internal energy is taken from the liberty internal_power tables, and leakage from the cell leakage attributes. The figures of merit are defined as`);
  EQ("T_{clk} = t_{comb} + t_{clk→Q} + t_{setup},     E_{op} = (E_{cycle} + P_{leak}T_{clk}) / ops", "11a");
  EQ("ADP = A·T_{clk},     EDP = E_{op}·T_{clk},     η = 1 / E_{op},     ρ = ops / (A·T_{clk})", "11b");
  P("where η is energy efficiency (TOPS/W, counting one MAC as one operation) and ρ is throughput density (GOPS/mm²). These are pre-layout results: wire capacitance and clock-tree power are excluded, which lowers absolute energy by an amount that is similar across architectures, so *relative* comparisons are more reliable than absolute values.");

  // ---------------- results: stand-alone ----------------
  H2("7.5 Stand-alone multiplier comparison (without HBF)");
  P("We first compare the multipliers as generic arithmetic units, independent of where they are deployed. Tables 8 and 9 report every figure of merit for the 8-bit designs in both libraries and both synthesis modes, and Fig. 8 normalises the main metrics to the array multiplier.");
  const SHORT = { bw_array8: "Array (BW)", bw_wallace8: "Wallace", bw_dadda8: "Dadda", bw_dadda8_ks: "Dadda + KS", booth8: "Booth-4", booth8_ks: "Booth-4 + KS",
    behav_mult8: "Tool a×b", bw_array8_p3: "Array, 3-stage", bw_dadda8_p2: "Dadda, 2-stage", booth8_p2: "Booth-4, 2-stage" };
  const fj = (pj) => g3(pj * 1000);
  const mrow = (r) => [SHORT[r.design] + (r.mode === "delay" ? " †" : ""), f(r.area, r.area < 100 ? 1 : 0), r.cells, f(r.tclk, 0), f(r.fmax_ghz, 2), fj(r.e_pj), fj(r.eu_pj), f(r.leak_uw, 3), g3(r.adp), g3(r.edp * 1000), f(r.gops_mm2, 0), g3(r.tops_w)];
  const MH = ["Design", "Area (µm²)", "Cells", "T_clk (ps)", "f_max (GHz)", "E/op Gauss. (fJ)", "E/op unif. (fJ)", "Leak (µW)", "ADP (µm²·ns)", "EDP (fJ·ns)", "GOPS /mm²", "TOPS /W"];
  const MW = [1580, 690, 520, 640, 640, 700, 700, 600, 760, 760, 760, 676];
  TBL("Signed 8×8 multipliers in SkyWater 130 nm (1.8 V). † = delay-driven synthesis (otherwise area-driven); KS = Kogge–Stone final adder. E/op includes leakage over one registered cycle; Gauss./unif. = operand statistics.",
    MH, T.mult8_sky130.map(mrow), MW, { size: 14 });
  TBL("Signed 8×8 multipliers in ASAP7 (0.7 V); notation as in Table 8.", MH, T.mult8_asap7.map(mrow), MW, { size: 14 });
  FIG("figures/fig8_mult_compare.png", "Six figures of merit for the 8-bit multiplier architectures in both libraries, normalised to the Baugh–Wooley array multiplier (area-driven synthesis, Gaussian operands). Values below 1 are better except for TOPS/W.");

  const bestE7 = best("e_pj", "asap7"), bestT7 = best("tclk", "asap7"), bestA7 = best("area", "asap7");
  const bestE1 = best("e_pj", "sky130"), bestT1 = best("tclk", "sky130"), bestA1 = best("area", "sky130");
  P(`**Area.** The smallest 8-bit multiplier is the ${bestA7.arch} in ASAP7 (${f(bestA7.area, 1)} µm²) and the ${bestA1.arch} in sky130 (${f(bestA1.area, 0)} µm²). The array multiplier occupies ${f(A.area, 1)} µm² (ASAP7) and ${f(As.area, 0)} µm² (sky130); Booth recoding ${Bo.area < A.area ? "reduces" : "increases"} area by ${f(Math.abs(Bo.area / A.area - 1) * 100, 0)}% in ASAP7 and ${Bs.area < As.area ? "reduces" : "increases"} it by ${f(Math.abs(Bs.area / As.area - 1) * 100, 0)}% in sky130: the saving of ${st.array8.fa - st.booth8.fa} full adders (Table 7) ${Bo.area < A.area ? "outweighs" : "is offset by"} the recoder and 5:1 selector logic. **Speed.** The fastest registered design is the ${bestT7.arch} (${f(bestT7.tclk, 0)} ps, ${f(bestT7.fmax_ghz, 2)} GHz) in ASAP7 and the ${bestT1.arch} (${f(bestT1.tclk, 0)} ps) in sky130, against ${f(A.tclk, 0)} ps and ${f(As.tclk, 0)} ps for the array. **Energy.** The lowest energy per multiply is achieved by the ${bestE7.arch} in ASAP7 (${f(bestE7.e_pj * 1000, 1)} fJ) and the ${bestE1.arch} in sky130 (${f(bestE1.e_pj, 3)} pJ); the array multiplier needs ${f(A.e_pj * 1000, 1)} fJ and ${f(As.e_pj, 3)} pJ respectively. The tool-default implementation of a×b (Yosys/ABC) is ${Be.area < A.area ? "smaller" : "larger"} than the hand-structured array (${f(Be.area, 1)} vs. ${f(A.area, 1)} µm²) and ${Be.tclk < A.tclk ? "faster" : "slower"} (${f(Be.tclk, 0)} vs. ${f(A.tclk, 0)} ps), which shows that a structural description still matters for arithmetic in an open-source flow.`);
  FIG("figures/fig9_width_scaling.png", "Scaling of area, critical path and energy with operand width (4, 8, 16 bits) in ASAP7.");
  const a16 = wid("Array (BW)", 16, "asap7"), d16 = wid("Dadda", 16, "asap7"), b16 = wid("Booth-4", 16, "asap7");
  const a4 = wid("Array (BW)", 4, "asap7"), a8w = wid("Array (BW)", 8, "asap7");
  const dk16 = T.width.find((r) => r.fam.startsWith("Dadda + Kogge") && r.lib === "asap7");
  const bk16 = T.width.find((r) => r.fam.startsWith("Booth-4 + Dadda + Kogge") && r.lib === "asap7");
  P(`**Width scaling (Fig. 9, Table 10).** Area and energy grow roughly quadratically with *n*: the ASAP7 array multiplier grows from ${f(a4.area, 1)} µm² at 4 bits to ${f(a8w.area, 1)} and ${f(a16.area, 1)} µm² at 8 and 16 bits. The critical path of the array grows linearly (${f(a4.delay, 0)}, ${f(a8w.delay, 0)}, ${f(a16.delay, 0)} ps), whereas tree multipliers with a parallel-prefix CPA grow logarithmically: at 16 bits Dadda+KS reaches ${f(dk16.delay, 0)} ps and Booth+KS ${f(bk16.delay, 0)} ps, i.e. ${f(a16.delay / Math.min(dk16.delay, bk16.delay), 2)}× faster than the array. At 16 bits Booth recoding also pays off in energy (${f(b16.e_pj * 1000, 1)} fJ vs. ${f(a16.e_pj * 1000, 1)} fJ for the array). For the 4- and 8-bit operands that dominate LLM inference, however, the array multiplier remains within ${f(worst8 * 100 - 100, 0)}% of the best design on area, clock period and energy in both libraries while being the most regular to lay out — a property that matters on a base die whose routing resources are shared with thousands of TSVs. Fig. 10 shows the complete energy–delay design space including both synthesis modes: delay-driven mapping shortens the clock period by ${f(dd.tmin, 0)}–${f(dd.tmax, 0)}% at up to ${f(dd.emax, 0)}% more energy per multiply, and the pipelined designs occupy a separate high-frequency, high-energy region.`);
  FIG("figures/fig10_pareto.png", "Energy–delay design space of all 8-bit multipliers (both synthesis modes; pipelined variants annotated). Points toward the lower left are Pareto-better.");
  const FAM = (x) => x.replace("Booth-4 + Dadda + Kogge-Stone CPA", "Booth-4 + KS").replace("Dadda + Kogge-Stone CPA", "Dadda + KS").replace("Array, 4-stage pipeline", "Array, 4-stage");
  const wrow = (r) => [FAM(r.fam), r.n, r.lib === "asap7" ? "ASAP7" : "sky130", f(r.area, r.area < 100 ? 1 : 0), f(r.delay, 0), f(r.fmax_ghz, 2), fj(r.e_pj), g3(r.edp * 1000), g3(r.tops_w)];
  TBL("Width sweep (4, 8, 16 bits) for all architectures in both libraries (area-driven synthesis, Gaussian operands). EDP in fJ·ns.",
    ["Architecture", "n", "Library", "Area (µm²)", "t_comb (ps)", "f_max (GHz)", "E/op (fJ)", "EDP (fJ·ns)", "TOPS/W"],
    T.width.map(wrow), [2300, 450, 850, 1000, 1000, 900, 900, 900, 726], { size: 16 });
  FIG("figures/fig11_operand_stats.png", "Energy per multiply for uniformly random and Gaussian (LLM-weight-like) operands.");
  const gu = (r) => r.e_pj / r.eu_pj;
  P(`**Operand statistics (Fig. 11).** LLM weights are concentrated near zero. In two's complement, small negative values have long runs of leading ones, so successive operands of opposite sign toggle the whole upper part of the partial-product array: with Gaussian operands the array multiplier consumes ${f(gu(A), 2)}× the energy it consumes for uniform operands in ASAP7 (${f(gu(As), 2)}× in sky130), whereas Booth recoding — whose recoder produces zero partial products for runs of identical bits — changes by ${f(gu(Bo), 2)}× (ASAP7). Data-aware choices such as sign-magnitude weight encoding or Booth recoding therefore deserve attention in an engine dedicated to LLM weights.`);
  FIG("figures/fig12_pipelining.png", "Effect of pipelining: clock-frequency gain, area overhead and energy per multiply relative to the combinational version.");
  const pp = (c, p, lib) => { const x = row(`mult8_${lib}`, c), y = row(`mult8_${lib}`, p); return [y.fmax_ghz / x.fmax_ghz, y.area / x.area, y.e_pj / x.e_pj]; };
  const [pf, pa, pe] = pp("bw_array8", "bw_array8_p3", "asap7");
  P(`**Pipelining (Fig. 12).** Three pipeline cuts raise the ASAP7 array multiplier's clock rate by ${f(pf, 2)}× at ${f((pa - 1) * 100, 0)}% more area and ${f((pe - 1) * 100, 0)}% ${pe >= 1 ? "more" : "less"} energy per multiply; the registers ${pe < 1 ? "also stop glitches from propagating through the array, which partly pays for their own switching energy" : "add clock and data switching energy that outweighs the glitch reduction"}. (Combinational designs are reported without input/output registers, whereas every pipeline register — including operand copies — is counted in the pipelined designs; liberty data put the clock-pin energy alone at ≈20 fJ per edge for a sky130 flip-flop.) Since a near-flash engine is sized by bandwidth rather than by clock rate (Section 7.6), pipelining is worthwhile only if it reduces the number of MAC units — and hence area — more than it adds per unit.`);
  const mr = (r) => [r.arch, r.lib === "asap7" ? "ASAP7" : "sky130", f(r.area, r.area < 100 ? 1 : 0), r.cells, f(r.tclk, 0), f(r.fmax_ghz, 2), fj(r.e_pj), f(r.leak_uw, 3), f(r.gops_mm2, 0), g3(r.tops_w)];
  TBL("MAC engines for the base die (area-driven synthesis, Gaussian operands; the systolic array performs 32 MACs per cycle and its E/op is per MAC).",
    ["Engine", "Library", "Area (µm²)", "Cells", "T_clk (ps)", "f_max (GHz)", "E/MAC (fJ)", "Leak (µW)", "GOPS/mm²", "TOPS/W"],
    T.macs.map(mr), [2200, 750, 800, 650, 750, 750, 800, 700, 850, 776], { size: 16 });
  const fp7 = mac("fp8_mac", "asap7"), pe7 = mac("int8_pe", "asap7"), sy7 = mac("systolic_4x8", "asap7");
  const fp1 = mac("fp8_mac", "sky130"), pe1 = mac("int8_pe", "sky130");
  P(`**MAC engines (Table 11).** The FP8 E4M3 MAC with exact 48-bit accumulation costs ${f(fp7.area / pe7.area, 2)}× the area and ${f(fp7.e_pj / pe7.e_pj, 2)}× the energy of the INT8 PE in ASAP7 (${f(fp1.area / pe1.area, 2)}× and ${f(fp1.e_pj / pe1.e_pj, 2)}× in sky130), although its significand multiplier is only 4×4: the alignment shifter and wide accumulator dominate. Its critical path is the stage that combines the barrel shifter with conditional negation (${f(fp1.tclk / 1000, 1)} ns in sky130), which is the obvious candidate for an extra pipeline cut. The 4×8 systolic array delivers ${f(sy7.gops_mm2, 0)} GOPS/mm² and ${g3(sy7.tops_w)} TOPS/W in ASAP7.`);

  // ---------------- results: system ----------------
  H2("7.6 System-level comparison: conventional HBM, HBF with far compute, and near-flash compute");
  P("We now place the engines in the three systems of Fig. 1 and Fig. 7: **S1**, a conventional package with weights in HBM and MACs on the accelerator die (no HBF); **S2**, HBF with weights streamed over the interposer to MACs on the accelerator; and **S3**, HBF with MACs in the base die. Table 12 lists the energy parameters. MAC energy, area and clock rate come from Tables 8–11 (accelerator MACs are taken from ASAP7; base-die MACs from either library, bracketing the unknown base-die process).");
  const D0 = T.defaults;
  TBL("System-model parameters (defaults; all are swept in Figs. 14–16).",
    ["Parameter", "Default", "Range explored", "Meaning"],
    [
      ["e_HBM", `${D0.e_hbm} pJ/bit`, "—", "HBM access + I/O to accelerator (HBM2-class [25])"],
      ["e_NAND", `${D0.e_nand} pJ/bit`, "0.5–10", "NAND sense + page-buffer read (assumption)"],
      ["e_TSV", `${D0.e_tsv} pJ/bit`, "—", "Core die → base die"],
      ["e_PHY", `${D0.e_phy} pJ/bit`, "0.25–3", "Base die → accelerator over interposer"],
      ["e_SRAM", `${D0.e_sram} pJ/bit`, "—", "Operand staging on the accelerator"],
      ["e_buf", `${D0.e_buf} pJ/bit`, "—", "Page buffer → near-flash MAC"],
      ["β_stack", `${D0.bw_stack / 1e12} TB/s`, "—", "External HBF read bandwidth per stack"],
      ["γ", `${D0.gamma}`, "1–4", "Internal-array / external bandwidth ratio"],
      ["P_budget", `${D0.p_budget} W`, "1–15", "Compute power allowed in one base die"],
      ["A_budget", `${D0.a_budget} mm²`, "—", "Base-die area available for compute"],
    ], [1400, 1500, 1500, 4626]);
  P("For a decode step with batch *B*, the weight-GEMV energy per generated token is");
  EQ("E_{token} = [ 8 W_{touched}(B) · e_{move}^{(S)} + B·P_{act}·E_{MAC} + E_{act}^{(S)} ] / B", "12a");
  EQ("e^{(S1)} = e_{HBM} + e_{SRAM},   e^{(S2)} = e_{NAND} + e_{TSV} + e_{PHY} + e_{SRAM},   e^{(S3)} = e_{NAND} + e_{TSV} + e_{buf}", "12b");
  P("where *E*_{act}^{(S3)} accounts for activations and FP32 partial sums that must still cross the interposer in S3. A base-die engine must sustain a MAC rate of γβ_{stack}·*B* (every weight byte is used once per sequence in the batch), which requires *U* = γβ_{stack}*B*/*f*_{clk} MAC units; the largest batch it can serve within the budgets is therefore");
  EQ("B_{max} = min( P_{budget} / (γ β_{stack} E_{MAC}),   A_{budget} · f_{clk} / (γ β_{stack} A_{MAC}) )", 13);
  FIG("figures/fig13_system_energy.png", "Weight-GEMV energy per generated token for S1 (HBM, no HBF), S2 (HBF, far compute) and S3 (HBF, near-flash MACs in a 130 nm or 7 nm base die); solid = data movement, hatched = MAC computation.");
  const e1 = sysE("MoE-1T", 1, "S1"), e2 = sysE("MoE-1T", 1, "S2"), e3a = sysE("MoE-1T", 1, "S3: HBF + near-flash MAC (7"), e3s = sysE("MoE-1T", 1, "S3: HBF + near-flash MAC (130");
  const e64_3s = sysE("MoE-1T", 64, "S3: HBF + near-flash MAC (130"), e64_2 = sysE("MoE-1T", 64, "S2"), e64_1 = sysE("MoE-1T", 64, "S1");
  P(`**Energy (Fig. 13, Table 13).** As in general-purpose computing [59], data movement dominates. For MoE-1T at *B* = 1, moving the touched weights costs ${f(e1.move, 1)} mJ per token in S1 and ${f(e2.move, 1)} mJ in S2, while the MACs themselves cost only ${f(e1.compute, 2)} mJ — a ratio of about ${f(e1.move / e1.compute, 0)}:1. S3 removes the interposer and accelerator-staging energy from the weight path and reduces the per-token energy to ${f(e3a.e_mJ, 1)} mJ with a 7 nm base die (${f((1 - e3a.e_mJ / e2.e_mJ) * 100, 0)}% below S2 and ${f((1 - e3a.e_mJ / e1.e_mJ) * 100, 0)}% below the conventional HBM system) and to ${f(e3s.e_mJ, 1)} mJ with a 130 nm base die. As the batch grows, weight traffic is amortised and MAC energy becomes visible: at *B* = 64 a 130 nm base-die engine costs ${f(e64_3s.compute, 2)} mJ per token in computation and the total is ${f(e64_3s.e_mJ, 2)} mJ versus ${f(e64_2.e_mJ, 2)} mJ for S2 and ${f(e64_1.e_mJ, 2)} mJ for S1, so **the multiplier's energy efficiency — and the base-die process — decides whether near-flash compute remains advantageous at large batch**.`);
  const er = (r) => [r.model, r.B, r.scn, f(r.move, r.move < 1 ? 3 : 2), f(r.compute, 3), f(r.act, 4), f(r.e_mJ, r.e_mJ < 1 ? 3 : 2)];
  TBL("Weight-GEMV energy per generated token (mJ) by system and batch.",
    ["Model", "B", "System", "Movement", "MAC", "Act./psum", "Total"], T.sys_energy.map(er), [1100, 500, 3800, 1000, 900, 900, 826], { size: 16 });
  FIG("figures/fig14_nmc_capacity.png", "(a) Largest decode batch a base-die engine can serve per stack (Equation (13)) for each multiplier architecture and library at 5 W and 20 mm²; (b) sensitivity of B_max to the thermal budget.");
  const c7 = cap("bw_array8", "asap7"), c1 = cap("bw_array8", "sky130"), cb7 = cap("booth8", "asap7"), cf7 = cap("fp8_mac", "asap7"), cf1 = cap("fp8_mac", "sky130");
  const bestCap = (lib) => T.nmc_capacity.filter((r) => r.lib === lib).reduce((b, r) => (r.Bmax > b.Bmax ? r : b));
  const bc7 = bestCap("asap7"), bc1 = bestCap("sky130");
  const lim = (r) => (r.Bmax_power < r.Bmax_area ? "power" : "area");
  P(`**How many MACs fit under the NAND? (Fig. 14, Table 14).** At 0.8 TB/s per stack every unit of batch requires ${f(c7.W_per_B, 3)} W and ${f(c7.mm2_per_B, 3)} mm² of INT8 array-multiplier PEs in ASAP7, but ${f(c1.W_per_B, 2)} W and ${f(c1.mm2_per_B, 2)} mm² in sky130. Under a 5 W, 20 mm² budget, a 7 nm base die can therefore serve decode batches up to *B*_{max} ≈ ${f(c7.Bmax, 0)} (limited by ${lim(c7)}), whereas a 130 nm base die is limited to *B*_{max} ≈ ${f(c1.Bmax, 1)} (limited by ${lim(c1)}). The best architecture for this purpose is the ${bc7.name} in ASAP7 (*B*_{max} ≈ ${f(bc7.Bmax, 0)}) and the ${bc1.name} in sky130 (*B*_{max} ≈ ${f(bc1.Bmax, 1)}). The FP8 MAC reduces *B*_{max} to ${f(cf7.Bmax, 0)} (ASAP7) and ${f(cf1.Bmax, 1)} (sky130). **The base-die process node, not the multiplier topology, is the first-order design decision**: it changes *B*_{max} by more than an order of magnitude, while the choice among multiplier architectures changes it by only ${f((Math.max(...T.nmc_capacity.filter(r => r.lib === "sky130" && r.design !== "fp8_mac").map(r => r.Bmax)) / Math.min(...T.nmc_capacity.filter(r => r.lib === "sky130" && r.design !== "fp8_mac").map(r => r.Bmax)) - 1) * 100, 0)}% (sky130) and ${f((Math.max(...T.nmc_capacity.filter(r => r.lib === "asap7" && r.design !== "fp8_mac").map(r => r.Bmax)) / Math.min(...T.nmc_capacity.filter(r => r.lib === "asap7" && r.design !== "fp8_mac").map(r => r.Bmax)) - 1) * 100, 0)}% (ASAP7), because the accumulator, pipeline and pass-through registers of a PE consume most of its energy.`);
  const cr = (r) => [r.name, r.lib === "asap7" ? "ASAP7" : "sky130", f(r.area, r.area < 100 ? 1 : 0), f(r.tclk, 0), f(r.e_mac_fJ, 1), f(r.W_per_B, 3), f(r.mm2_per_B, 3), f(r.Bmax_power, 1), f(r.Bmax_area, 1), f(r.Bmax, 1)];
  TBL("Near-flash MAC capacity per HBF stack (0.8 TB/s, γ = 1, 5 W, 20 mm²). PE = multiplier + shared accumulator/registers.",
    ["PE multiplier", "Library", "PE area (µm²)", "T_clk (ps)", "E/MAC (fJ)", "W per unit B", "mm² per unit B", "B_max (power)", "B_max (area)", "B_max"],
    T.nmc_capacity.map(cr), [1700, 800, 900, 800, 800, 900, 900, 800, 800, 626], { size: 16 });
  FIG("figures/fig15_sensitivity.png", "Sensitivity of near-flash compute (S3) energy per token relative to S2 (left) and to the conventional HBM system S1 (right), over NAND array read energy and interposer PHY energy (MoE-1T, B = 1, 7 nm base die). Values below 1 favour near-flash compute.");
  const H = T.heat, i2 = H.e_nand.indexOf(2), j05 = H.e_phy.indexOf(0.5);
  const maxS1 = Math.max(...H.S1.flat()), minS1 = Math.min(...H.S1.flat());
  const crossNand = H.e_nand.filter((e, i) => H.S1[i][j05] < 1);
  P(`**Sensitivity (Fig. 15).** Relative to S2, near-flash compute is always beneficial — the ratio ranges from ${f(Math.min(...H.S2.flat()), 2)} to ${f(Math.max(...H.S2.flat()), 2)} over the explored grid — because it removes PHY and staging energy without adding any. Relative to the conventional HBM system the answer depends on the NAND read energy: at the default 0.5 pJ/bit PHY, S3 is better than S1 whenever *e*_{NAND} ≤ ${f(Math.max(...crossNand), 1)} pJ/bit, and the ratio spans ${f(minS1, 2)}–${f(maxS1, 2)} over the grid. The energy-per-bit of the NAND read path is thus the single most important unknown for HBF as an energy-efficient — not merely capacity-efficient — technology.`);
  FIG("figures/fig16_nmc_throughput.png", "MoE-1T decode throughput per package (TPOT ≤ 50 ms, 32k context) when near-flash MAC engines consume weights at γ times the external HBF bandwidth, limited to the batch B_max they can serve. Dotted line: conventional HBM-8 system on four packages.");
  const th = (cfg, N, lib, g) => T.nmc_thr.find((r) => r.cfg === cfg && r.N === N && r.lib === lib && r.gamma === g);
  const t41 = th("HBM-4/HBF-4", 4, "asap7", 1), t44 = th("HBM-4/HBF-4", 4, "asap7", 4), t11 = th("HBM-4/HBF-4", 1, "asap7", 1), t14 = th("HBM-4/HBF-4", 1, "asap7", 4);
  const ts44 = th("HBM-4/HBF-4", 4, "sky130", 4), t62 = th("HBM-2/HBF-6", 4, "asap7", 2);
  P(`**Throughput (Fig. 16).** If the base die can read the arrays faster than its external interface allows (γ > 1) — plausible because internal TSV bandwidth is not limited by interposer beachfront — near-flash engines convert that internal bandwidth directly into decode throughput. For MoE-1T on four HBM-4/HBF-4 packages, throughput rises from ${f(t41.thr, 0)} tokens/s per package at γ = 1 to ${f(t44.thr, 0)} at γ = 4 with 7 nm engines, compared with ${f(T.nmc_thr_base["4"], 0)} for the conventional HBM-8 system on the same four packages; a single package reaches ${f(t14.thr, 0)} tokens/s at γ = 4 (from ${f(t11.thr, 0)}). With 130 nm engines the batch limit of Equation (13) caps the gain at ${f(ts44.thr, 0)} tokens/s per package even at γ = 4. The curves are not always monotonic in γ: each unit of batch costs γ times more MAC power (Equation (13)), so a faster array read shrinks the batch the base die can serve. Because the accelerator's own MACs remain available, a real system would fall back to far compute (S2) for batches above *B*_{max}; the near-flash path should be read as an additional, energy-efficient mode for small batches rather than a replacement. Near-flash multipliers are therefore an *enabler* of the HBF capacity advantage only when the base die is fabricated in a node advanced enough to host thousands of MAC units within a few watts.`);

  H2("7.7 Design guidelines");
  NUM([
    "**Choose the base-die node first.** It changes the sustainable near-flash batch by more than an order of magnitude (Table 14); multiplier topology changes it by at most about 20%.",
    `**Use array or Booth multipliers at INT8/FP8 widths.** At ≤ 8 bits the regular array multiplier is within ${f(worst8 * 100 - 100, 0)}% of the best design on area, clock period and energy and is the easiest to place between TSV fields; Booth recoding (with a parallel-prefix CPA) is the best single choice in ASAP7 and becomes clearly superior at 16 bits.`,
    "**Optimise the accumulator and registers, not only the multiplier.** In a PE they dominate energy; narrower accumulators, clock gating and fewer pass-through registers yield more than a better multiplier.",
    "**Add a parallel-prefix CPA where it shrinks the engine.** Kogge–Stone shortens the critical path; in a bandwidth-sized engine a faster clock means fewer MAC units, so the CPA pays off only if the area it saves per unit of batch exceeds its own cost.",
    "**Pipeline FP8 alignment, not the multiplier.** In the FP8 MAC the shifter-plus-negation stage, not the 4×4 significand multiplier, limits the clock.",
    "**Keep only GEMV in the base die.** Attention over the KV cache stays in HBM next to the accelerator (Table 1); the near-flash engine serves dense and expert FFN weights, where the weight-to-activation ratio is largest.",
  ], "guide");
};
