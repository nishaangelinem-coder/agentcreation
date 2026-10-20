#!/usr/bin/env python3
"""Download the open standard-cell libraries, merge ASAP7, and generate Verilog cell models."""
import os, re, subprocess, urllib.request, gzip, shutil
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("lib", exist_ok=True); os.makedirs("syn/cells", exist_ok=True)
ORFS = "https://raw.githubusercontent.com/The-OpenROAD-Project/OpenROAD-flow-scripts/master/flow/platforms"
FILES = {
    "sky130_fd_sc_hd__tt_025C_1v80.lib": f"{ORFS}/sky130hd/lib/sky130_fd_sc_hd__tt_025C_1v80.lib",
    "asap7sc7p5t_SIMPLE_RVT_TT_nldm_211120.lib": f"{ORFS}/asap7/lib/NLDM/asap7sc7p5t_SIMPLE_RVT_TT_nldm_211120.lib.gz",
    "asap7sc7p5t_AO_RVT_TT_nldm_211120.lib": f"{ORFS}/asap7/lib/NLDM/asap7sc7p5t_AO_RVT_TT_nldm_211120.lib.gz",
    "asap7sc7p5t_OA_RVT_TT_nldm_211120.lib": f"{ORFS}/asap7/lib/NLDM/asap7sc7p5t_OA_RVT_TT_nldm_211120.lib.gz",
    "asap7sc7p5t_INVBUF_RVT_TT_nldm_220122.lib": f"{ORFS}/asap7/lib/NLDM/asap7sc7p5t_INVBUF_RVT_TT_nldm_220122.lib.gz",
    "asap7sc7p5t_SEQ_RVT_TT_nldm_220123.lib": f"{ORFS}/asap7/lib/NLDM/asap7sc7p5t_SEQ_RVT_TT_nldm_220123.lib",
}
for name, url in FILES.items():
    dst = f"lib/{name}"
    if os.path.exists(dst): continue
    tmp = dst + (".gz" if url.endswith(".gz") else "")
    urllib.request.urlretrieve(url, tmp)
    if url.endswith(".gz"):
        with gzip.open(tmp, "rb") as fi, open(dst, "wb") as fo: shutil.copyfileobj(fi, fo)
        os.remove(tmp)

def split(txt):
    i = re.search(r"\n\s*cell\s*\(", txt).start(); j = txt.rstrip().rfind("}")
    return txt[:i], txt[i:j]
asap = [n for n in FILES if n.startswith("asap7")]
head0, _ = split(open(f"lib/{asap[0]}").read())
templates, cells = {}, []
for n in asap:
    h, c = split(open(f"lib/{n}").read())
    for m in re.finditer(r"\n(\s*lu_table_template\s*\((\w+)\)\s*\{.*?\})", h, re.S):
        templates.setdefault(m.group(2), m.group(1))
    cells.append(c)
head = re.sub(r"\n\s*lu_table_template\s*\(\w+\)\s*\{.*?\}", "", head0, flags=re.S)
head = head.replace("library (asap7sc7p5t_SIMPLE_RVT_TT_nldm_211120)", "library (asap7_rvt_tt_merged)")
open("lib/asap7_rvt_tt.lib", "w").write(head + "\n" + "\n".join(templates.values()) + "\n" + "\n".join(cells) + "\n}\n")
for lib, out in (("lib/asap7_rvt_tt.lib", "syn/cells/asap7_cells.v"), ("lib/sky130_fd_sc_hd__tt_025C_1v80.lib", "syn/cells/sky130_cells.v")):
    subprocess.run(["yosys", "-q", "-p", f"read_liberty -ignore_miss_func {lib}; write_verilog -noattr {out}"], check=True)
print("libraries ready")
