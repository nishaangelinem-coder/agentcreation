#!/usr/bin/env python3
"""Synthesis + glitch-aware power flow for the multiplier / MAC study.

For every design x library x optimisation goal:
  1. Yosys: synth -flatten, dfflibmap, ABC liberty mapping (area- or delay-driven), stat.
  2. Critical path: ABC `stime -p` (+ flop clk->Q and setup from the liberty for sequential designs).
  3. Power: unit-delay gate-level simulation (Icarus) of the mapped netlist with two operand
     statistics (uniform random, and Gaussian 'LLM-weight-like'), VCD toggle counting,
     switching energy 0.5*C*V^2 from liberty pin capacitances, internal energy from
     liberty internal_power tables (median entry), leakage from liberty.
Results are written to results/syn_results.json. Pre-layout: no wire capacitance, ideal clock.
"""
import json, os, re, subprocess, sys, random, statistics, math, functools

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
B = "build/syn"
os.makedirs(B, exist_ok=True)

LIBS = {
    "sky130": dict(lib="lib/sky130_fd_sc_hd__tt_025C_1v80.lib", cells="syn/cells/sky130_cells.v",
                   V=1.8, cap_unit=1e-12, energy_unit=1e-12, leak_unit=1e-9, time_unit=1e-9, node="SkyWater 130 nm"),
    "asap7": dict(lib="lib/asap7_rvt_tt.lib", cells="syn/cells/asap7_cells.v",
                  V=0.7, cap_unit=1e-15, energy_unit=1e-15, leak_unit=1e-12, time_unit=1e-12, node="ASAP7 7 nm"),
}

def mult(name, n, lat=0, arch=None):
    return dict(top=name, files=[f"rtl/{name}.v"], kind="mult", n=n, lat=lat, ops=1, arch=arch)

DESIGNS = {}
for n in (4, 8, 16):
    DESIGNS[f"bw_array{n}"] = mult(f"bw_array{n}", n, arch="Array (Baugh-Wooley)")
    DESIGNS[f"bw_wallace{n}"] = mult(f"bw_wallace{n}", n, arch="Wallace tree")
    DESIGNS[f"bw_dadda{n}"] = mult(f"bw_dadda{n}", n, arch="Dadda tree")
    DESIGNS[f"booth{n}"] = mult(f"booth{n}", n, arch="Radix-4 Booth + Dadda")
    DESIGNS[f"behav_mult{n}"] = mult(f"behav_mult{n}", n, arch="Tool default (a*b)")
for n in (8, 16):
    DESIGNS[f"bw_dadda{n}_ks"] = mult(f"bw_dadda{n}_ks", n, arch="Dadda + Kogge-Stone CPA")
    DESIGNS[f"booth{n}_ks"] = mult(f"booth{n}_ks", n, arch="Booth-4 + Dadda + Kogge-Stone CPA")
DESIGNS["bw_array8_p3"] = mult("bw_array8_p3", 8, 3, "Array, 3-stage pipeline")
DESIGNS["bw_dadda8_p2"] = mult("bw_dadda8_p2", 8, 2, "Dadda, 2-stage pipeline")
DESIGNS["booth8_p2"] = mult("booth8_p2", 8, 2, "Booth-4, 2-stage pipeline")
DESIGNS["bw_array16_p4"] = mult("bw_array16_p4", 16, 4, "Array, 4-stage pipeline")
DESIGNS["u_array4"] = dict(top="u_array4", files=["rtl/u_array4.v"], kind="umult", n=4, lat=0, ops=1, arch="Unsigned 4x4 array (FP8 mantissa)")
DESIGNS["int8_pe"] = dict(top="int8_pe", files=["rtl/int8_pe.v", "rtl/bw_array8.v"], kind="pe", n=8, lat=2, ops=1, arch="INT8 MAC PE (array mult + 32b acc)")
DESIGNS["fp8_mac"] = dict(top="fp8_mac", files=["rtl/fp8_mac.v", "rtl/u_array4.v"], kind="fp8", n=8, lat=3, ops=1, arch="FP8 E4M3 MAC (exact 48b acc)")
DESIGNS["systolic_4x8"] = dict(top="systolic_array", files=["rtl/systolic_array.v", "rtl/int8_pe.v", "rtl/bw_array8.v"],
                              kind="sys", n=8, lat=2, ops=32, arch="4x8 INT8 systolic array")

CYCLES = int(os.environ.get("SYN_CYCLES", 3000))
PERIOD_PS = 20000  # simulation period; unit (1 ps) cell delays settle well within it

# ---------------------------------------------------------------- liberty parsing
@functools.lru_cache(None)
def lib_text(path):
    return open(path).read()

@functools.lru_cache(None)
def cell_index(path):
    t = lib_text(path)
    return {m.group(1): m.start() for m in re.finditer(r'\bcell\s*\(\s*"?([\w]+)"?\s*\)\s*\{', t)}

def block(t, i):
    j = t.index("{", i); d = 0
    for k in range(j, len(t)):
        if t[k] == "{": d += 1
        elif t[k] == "}":
            d -= 1
            if d == 0: return t[j + 1:k]
    raise ValueError

def sub_blocks(t, kw):
    out = []
    for m in re.finditer(rf'\b{kw}\s*\(\s*"?([\w]*)"?\s*\)\s*\{{', t):
        out.append((m.group(1), block(t, m.start())))
    # keep only top-level ones (not nested inside a previous one)
    res, end = [], -1
    for m, (nm, b) in zip(re.finditer(rf'\b{kw}\s*\(\s*"?([\w]*)"?\s*\)\s*\{{', t), out):
        if m.start() > end:
            res.append((nm, b)); end = m.start() + len(b)
    return res

def nums(s):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s)]

def table_median(b, kw):
    vals = []
    for _, tb in sub_blocks(b, kw):
        m = re.search(r"values\s*\((.*?)\)\s*;", tb, re.S)
        if m: vals += nums(m.group(1))
    return statistics.median(vals) if vals else None

@functools.lru_cache(None)
def cell_info(path, cell):
    t = lib_text(path)
    cb = block(t, cell_index(path)[cell])
    area = float(re.search(r"\barea\s*:\s*([\d.eE+-]+)", cb).group(1))
    m = re.search(r"cell_leakage_power\s*:\s*([\d.eE+-]+)", cb)
    if m:
        leak = float(m.group(1))
    else:
        lv = [float(x) for x in re.findall(r"leakage_power\s*\(\s*\)\s*\{[^}]*?value\s*:\s*\"?([\d.eE+-]+)", cb)]
        leak = statistics.mean(lv) if lv else 0.0
    pins = {}
    for pn, pb in sub_blocks(cb, "pin"):
        d = re.search(r"direction\s*:\s*\"?(\w+)", pb)
        c = re.search(r"\bcapacitance\s*:\s*([\d.eE+-]+)", pb)
        eint = []
        for _, ib in sub_blocks(pb, "internal_power"):
            r, f = table_median(ib, "rise_power"), table_median(ib, "fall_power")
            v = [x for x in (r, f) if x is not None]
            if v: eint.append(sum(v) / len(v))
        clkq = setup = None
        for _, tb in sub_blocks(pb, "timing"):
            rp = re.search(r'related_pin\s*:\s*"?(\w+)', tb)
            tt = re.search(r'timing_type\s*:\s*"?(\w+)', tb)
            if tt and tt.group(1) == "rising_edge":
                v = [x for x in (table_median(tb, "cell_rise"), table_median(tb, "cell_fall")) if x is not None]
                clkq = max(v) if v else None
            if tt and tt.group(1) == "setup_rising":
                v = [x for x in (table_median(tb, "rise_constraint"), table_median(tb, "fall_constraint")) if x is not None]
                setup = max(v) if v else None
        pins[pn] = dict(dir=d.group(1) if d else "input", cap=float(c.group(1)) if c else 0.0,
                        eint=max(0.0, statistics.mean(eint)) if eint else 0.0, clkq=clkq, setup=setup)
    return dict(area=area, leak=leak, pins=pins)

# ---------------------------------------------------------------- stimulus
def e4m3_table():
    vals = []
    for v in range(256):
        s, e, m = v >> 7, (v >> 3) & 15, v & 7
        x = (m / 8) * 2 ** -6 if e == 0 else (1 + m / 8) * 2 ** (e - 7)
        vals.append(-x if s else x)
    return vals
E4M3 = e4m3_table()

def to_e4m3(x):
    return min(range(256), key=lambda i: abs(E4M3[i] - x))

def operand(n, profile, rng):
    if profile == "uniform":
        return rng.randrange(1 << n)
    sigma = (1 << (n - 1)) / 8.0          # LLM-like: concentrated near zero, ~3% of full scale
    v = max(-(1 << (n - 1)), min((1 << (n - 1)) - 1, int(round(rng.gauss(0, sigma)))))
    return v & ((1 << n) - 1)

def fp8_operand(profile, rng):
    return rng.randrange(256) if profile == "uniform" else to_e4m3(rng.gauss(0, 1.0))

def stimulus(d, profile, seed=1):
    """Returns (list of (port, width)) and per-cycle values."""
    rng = random.Random(seed)
    k, n = d["kind"], d["n"]
    if k in ("mult", "umult"):
        prof = "uniform" if k == "umult" else profile
        ports = [("a", n), ("b", n)]
        vec = [[operand(n, prof, rng), operand(n, prof, rng)] for _ in range(CYCLES)]
    elif k == "pe":
        ports = [("clr", 1), ("a_in", 8), ("w_in", 8)]
        vec = [[int(c % 64 == 0), operand(8, profile, rng), operand(8, profile, rng)] for c in range(CYCLES)]
    elif k == "fp8":
        ports = [("rst", 1), ("clr", 1), ("valid", 1), ("w", 8), ("x", 8)]
        vec = [[int(c < 2), int(c % 64 == 2), 1, fp8_operand(profile, rng), fp8_operand(profile, rng)] for c in range(CYCLES)]
    elif k == "sys":
        ports = [("clr", 1), ("a_west", 32), ("w_north", 64)]
        vec = []
        for c in range(CYCLES):
            a = sum(operand(8, profile, rng) << (8 * i) for i in range(4))
            w = sum(operand(8, profile, rng) << (8 * i) for i in range(8))
            vec.append([int(c < 16), a, w])
    return ports, vec

# ---------------------------------------------------------------- flow steps
def run(cmd, log=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if log:
        open(log, "w").write(r.stdout + r.stderr)
    if r.returncode:
        sys.stderr.write(r.stdout[-3000:] + r.stderr[-3000:])
        raise RuntimeError(cmd)
    return r.stdout

def synth(dname, d, lname, mode):
    L = LIBS[lname]; tag = f"{dname}.{lname}.{mode}"
    D = "-D 1" if mode == "delay" else ""
    script = (f"+strash;&get,-n;&fraig,-x;&put;dc2;strash;&get,-n;&dch,-f;&nf,{{D}};&put;"
              + ("buffer;upsize,{D};dnsize,{D};" if mode == "delay" else "") + "stime,-p")
    ys = f"""read_verilog {' '.join(d['files'])}
synth -flatten -top {d['top']}
dfflibmap -liberty {L['lib']}
abc -liberty {L['lib']} {D} -script "{script}"
opt_clean -purge
stat -liberty {L['lib']}
write_verilog -noattr {B}/{tag}.v
write_json {B}/{tag}.json
"""
    open(f"{B}/{tag}.ys", "w").write(ys)
    log = run(f"yosys -q -l {B}/{tag}.log -s {B}/{tag}.ys >/dev/null; cat {B}/{tag}.log")
    area = float(re.findall(r"Chip area for (?:top )?module.*?:\s*([\d.]+)", log)[-1])
    ncell = int(re.findall(r"Number of cells:\s*(\d+)", log)[-1])
    delays = [float(x) for x in re.findall(r"Delay\s*=\s*([\d.]+)\s*ps", log)]
    delay_ps = max(delays) if delays else 0.0
    return tag, area, ncell, delay_ps

def unit_delay_models(lname):
    out = f"{B}/{lname}_cells_ud.v"
    if os.path.exists(out): return out
    src = open(LIBS[lname]["cells"]).read()
    res = []
    for mod in re.split(r"(?=\nmodule )", src):
        outs = re.findall(r"^\s*output\s+(?:\[[^\]]*\]\s*)?(\w+)\s*;", mod, re.M)
        for o in outs:
            mod = re.sub(rf"^(\s*)assign\s+{o}\s*=", rf"\1assign #1 {o} =", mod, flags=re.M)
        res.append(mod)
    open(out, "w").write("`timescale 1ps/1ps\n" + "".join(res))
    return out

def simulate(tag, d, lname, profile):
    ports, vec = stimulus(d, profile)
    seq = d["kind"] != "mult" or d["lat"] > 0
    nl = json.load(open(f"{B}/{tag}.json"))["modules"][d["top"]]
    has_clk = "clk" in nl["ports"]
    width = sum(w for _, w in ports)
    stim = f"{B}/{tag}.{profile}.hex"
    with open(stim, "w") as f:
        for v in vec:
            word = 0
            for (p, w), x in zip(ports, v):
                word = (word << w) | (x & ((1 << w) - 1))
            f.write(f"{word:x}\n")
    conns, off = [], width
    for p, w in ports:
        off -= w
        conns.append(f".{p}(stim[{off + w - 1}:{off}])")
    outs = [p for p, info in nl["ports"].items() if info["direction"] == "output"]
    owidth = {p: len(nl["ports"][p]["bits"]) for p in outs}
    for p in outs:
        conns.append(f".{p}(o_{p})")
    tb = f"""`timescale 1ps/1ps
module tbp;
  reg clk = 0; reg [{width-1}:0] mem [0:{CYCLES-1}]; reg [{width-1}:0] stim = 0; integer c;
  {' '.join(f'wire [{owidth[p]-1}:0] o_{p};' for p in outs)}
  {d['top']} dut ({'.clk(clk), ' if has_clk else ''}{', '.join(conns)});
  initial begin
    $readmemh("{stim}", mem);
    $dumpfile("{B}/{tag}.{profile}.vcd"); $dumpvars(1, dut);
    for (c = 0; c < {CYCLES}; c = c + 1) begin
      stim = mem[c]; #{PERIOD_PS // 2} clk = 1; #{PERIOD_PS // 2} clk = 0;
    end
    $finish;
  end
endmodule
"""
    open(f"{B}/{tag}.tbp.v", "w").write(tb)
    run(f"iverilog -g2012 -o {B}/{tag}.simp {B}/{tag}.tbp.v {B}/{tag}.v {unit_delay_models(lname)} && vvp -n {B}/{tag}.simp >/dev/null")
    return toggle_counts(f"{B}/{tag}.{profile}.vcd", nl)

def toggle_counts(vcd, nl):
    """Per netlist bit-id toggle counts from a VCD dumped at the DUT scope."""
    code2bits, last, tog = {}, {}, {}
    name2bits = {n: info["bits"] for n, info in nl["netnames"].items()}
    with open(vcd) as f:
        defs = True
        for line in f:
            if defs:
                if line.startswith("$var"):
                    p = line.split()
                    code, nm = p[3], p[4]
                    if nm in name2bits:
                        code2bits[code] = name2bits[nm]
                elif line.startswith("$enddefinitions"):
                    defs = False
                continue
            if not line or line[0] in "#$": continue
            if line[0] in "01xz":
                val, code = line[0], line[1:].strip()
                vals = [val]
            elif line[0] == "b":
                v, code = line[1:].split()
                vals = None
                bits = code2bits.get(code)
                if bits is None: continue
                v = v.rjust(len(bits), v[0] if v[0] in "xz" else "0")
                vals = list(reversed(v))     # LSB first, matching JSON bit order
            else:
                continue
            bits = code2bits.get(code)
            if bits is None: continue
            prev = last.get(code)
            for i, (b, nv) in enumerate(zip(bits, vals)):
                if prev is not None and prev[i] != nv and prev[i] in "01" and nv in "01" and isinstance(b, int):
                    tog[b] = tog.get(b, 0) + 1
            last[code] = vals
    return tog

def power(tag, d, lname, tog):
    L = LIBS[lname]; lib = L["lib"]
    nl = json.load(open(f"{B}/{tag}.json"))["modules"][d["top"]]
    load = {}; e_int = 0.0; leak = 0.0; ffs = []
    dff_d_cap = None
    for cname, c in nl["cells"].items():
        ci = cell_info(lib, c["type"])
        leak += ci["leak"]
        for pin, bits in c["connections"].items():
            pi = ci["pins"].get(pin)
            if not pi: continue
            for b in bits:
                if not isinstance(b, int): continue
                if pi["dir"] == "input":
                    load[b] = load.get(b, 0.0) + pi["cap"]
                e_int += tog.get(b, 0) * pi["eint"]
            if pi["clkq"] is not None:
                ffs.append(ci)
            if pi.get("setup") is not None and dff_d_cap is None:
                dff_d_cap = pi["cap"]
    for p, info in nl["ports"].items():   # primary outputs drive one flop D pin
        if info["direction"] == "output":
            for b in info["bits"]:
                if isinstance(b, int):
                    load[b] = load.get(b, 0.0) + (dff_d_cap or 0.0)
    e_sw = sum(t * 0.5 * load.get(b, 0.0) * L["V"] ** 2 for b, t in tog.items())
    e_cycle = (e_sw * L["cap_unit"] + e_int * L["energy_unit"]) / CYCLES   # joules
    return e_cycle, leak * L["leak_unit"], ffs

def seq_overhead(ffs, lname):
    if not ffs: return 0.0
    ci = ffs[0]
    clkq = max(p["clkq"] for p in ci["pins"].values() if p["clkq"] is not None)
    setup = max([p["setup"] for p in ci["pins"].values() if p["setup"] is not None] or [0.0])
    scale = LIBS[lname]["time_unit"] / 1e-12
    return (clkq + setup) * scale

def main(only=None):
    out = {}
    path = os.environ.get("SYN_OUT", "results/syn_results.json")
    if os.path.exists(path): out = json.load(open(path))
    for dname, d in DESIGNS.items():
        if only and dname not in only: continue
        for lname in [l for l in LIBS if l in os.environ.get("SYN_LIBS", "sky130,asap7").split(",")]:
            for mode in ("area", "delay"):
                tag, area, ncell, dly = synth(dname, d, lname, mode)
                rec = dict(design=dname, arch=d["arch"], n=d["n"], lat=d["lat"], ops=d["ops"], lib=lname, mode=mode,
                           area_um2=area, cells=ncell, comb_delay_ps=dly)
                for prof in ("uniform", "gaussian"):
                    tog = simulate(tag, d, lname, prof)
                    e_cycle, p_leak, ffs = power(tag, d, lname, tog)
                    rec[f"e_cycle_J_{prof}"] = e_cycle
                    rec["p_leak_W"] = p_leak
                    os.remove(f"{B}/{tag}.{prof}.vcd")
                ovh = seq_overhead(ffs, lname) if (d["lat"] or d["kind"] != "mult") else 0.0
                rec["seq_overhead_ps"] = ovh
                rec["t_clk_ps"] = dly + ovh
                out[tag] = rec
                print(f"{tag:32s} area={area:10.2f} cells={ncell:6d} delay={dly:8.1f}ps "
                      f"E/cyc(u)={rec['e_cycle_J_uniform']*1e12:8.4f}pJ E/cyc(g)={rec['e_cycle_J_gaussian']*1e12:8.4f}pJ "
                      f"leak={p_leak*1e6:8.3f}uW", flush=True)
                json.dump(out, open(path, "w"), indent=1)

if __name__ == "__main__":
    main(sys.argv[1:] or None)
