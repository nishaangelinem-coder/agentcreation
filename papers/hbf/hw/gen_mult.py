#!/usr/bin/env python3
"""Generate structural Verilog multipliers.

  arch   : 'array' (carry-save array), 'wallace', 'dadda' (Baugh-Wooley partial products)
           or 'booth' (radix-4 Booth partial products + Dadda tree)
  signed : True -> modified Baugh-Wooley two's-complement; False -> unsigned
  pipe   : list of reduction-step indices after which a register stage is inserted
           (array: after CSA row k; dadda: after tree level k); 'F' splits the final adder

Usage: gen_mult.py NAME N {array|dadda} {s|u} [stage,...]
"""
import sys

class Gen:
    def __init__(self, name, n, arch, signed, pipe):
        arch, _, cpa = arch.partition("+")
        self.cpa = cpa or "rca"
        self.name, self.n, self.arch, self.signed = name, n, arch, signed
        self.pipe = set(pipe)
        self.W = 2 * n
        self.lines, self.cnt, self.regs = [], 0, []
        self.stage = 0
        self.depth = {}               # signal -> logic depth in 2-input-gate levels
        self.nfa = self.nha = 0
        self.stage_depths = []
        self.out_bits = []

    def w(self, expr, cost=1, deps=()):
        self.cnt += 1
        nm = f"n{self.cnt}"
        self.lines.append(f"  wire {nm} = {expr};")
        self.depth[nm] = cost + max([self.depth.get(x, 0) for x in deps] or [0])
        return nm

    def fa(self, a, b, c):
        self.nfa += 1
        s = self.w(f"{a} ^ {b} ^ {c}", 2, (a, b, c))
        co = self.w(f"({a} & {b}) | ({c} & ({a} ^ {b}))", 2, (a, b, c))
        return s, co

    def ha(self, a, b):
        self.nha += 1
        return self.w(f"{a} ^ {b}", 1, (a, b)), self.w(f"{a} & {b}", 1, (a, b))

    def reg(self, sig):
        """Insert a pipeline register for one signal (constants pass through)."""
        if sig in ("1'b0", "1'b1"):
            return sig
        self.cnt += 1
        nm = f"r{self.cnt}"
        self.regs.append((nm, sig))
        self.stage_depths.append(self.depth.get(sig, 0))
        self.depth[nm] = 0
        return nm

    def reg_ops(self, a, b):
        return [self.reg(x) for x in a], [self.reg(x) for x in b]

    def pp(self, a, b, i, j):
        n = self.n
        e = f"{a[j]} & {b[i]}"
        if self.signed and ((i == n - 1) != (j == n - 1)):
            return self.w(f"~({e})", 1, (a[j], b[i]))
        return self.w(e, 1, (a[j], b[i]))

    def pp_row(self, a, b, i):
        row = {i + j: self.pp(a, b, i, j) for j in range(self.n)}
        if self.signed and i == 0:
            row[self.n] = "1'b1"
            row[self.W - 1] = "1'b1"
        return row

    # ---------------- reductions ----------------
    def array(self, a, b):
        S, C = self.pp_row(a, b, 0), self.pp_row(a, b, 1)
        for k in range(2, self.n):
            row = self.pp_row(a, b, k)
            nS, nC = {}, {}
            for c in range(self.W):
                bits = [x for x in (S.get(c), C.get(c), row.get(c)) if x]
                if len(bits) == 3:
                    s, co = self.fa(*bits)
                elif len(bits) == 2:
                    s, co = self.ha(*bits)
                elif len(bits) == 1:
                    s, co = bits[0], None
                else:
                    continue
                nS[c] = s
                if co and c + 1 < self.W:
                    nC[c + 1] = co
            S, C = nS, nC
            if k in self.pipe:
                S = {c: self.reg(v) for c, v in S.items()}
                C = {c: self.reg(v) for c, v in C.items()}
                a, b = self.reg_ops(a, b)
        return S, C

    def bw_cols(self, a, b):
        cols = [[] for _ in range(self.W)]
        for i in range(self.n):
            for c, v in self.pp_row(a, b, i).items():
                cols[c].append(v)
        return cols

    def booth_cols(self, a, b):
        """Radix-4 Booth recoding of b (signed); n/2 rows of n+1 bits.
        Row value = ppbits (n+1-bit two's complement) + neg; the sign bit s is
        replaced by ~s at weight 2^(2i+n) plus a global constant -sum 2^(2i+n)."""
        n, W = self.n, self.W
        assert self.signed and n % 2 == 0
        cols = [[] for _ in range(W)]
        ax = a + [a[n - 1]]                       # sign-extend a to n+1 bits
        const = 0
        for i in range(n // 2):
            b_hi, b_mid = b[2 * i + 1], b[2 * i]
            b_lo = b[2 * i - 1] if i > 0 else "1'b0"
            one = self.w(f"{b_mid} ^ {b_lo}", 1, (b_mid, b_lo))
            two = self.w(f"({b_hi} & ~{b_mid} & ~{b_lo}) | (~{b_hi} & {b_mid} & {b_lo})", 2, (b_hi, b_mid, b_lo))
            neg = b_hi
            for k in range(n + 1):
                prev = ax[k - 1] if k > 0 else "1'b0"
                bit = self.w(f"(({one} & {ax[k]}) | ({two} & {prev})) ^ {neg}", 3, (one, two, ax[k], prev, neg))
                col = 2 * i + k
                if k == n:
                    bit = self.w(f"~{bit}", 0, (bit,))
                    const -= 1 << col
                if col < W:
                    cols[col].append(bit)
            cols[2 * i].append(neg)
        const %= 1 << W
        for c in range(W):
            if (const >> c) & 1:
                cols[c].append("1'b1")
        return cols

    def wallace(self, a, b):
        cols = self.booth_cols(a, b) if self.arch == "booth" else self.bw_cols(a, b)
        lvl = 0
        while max(len(c) for c in cols) > 2:
            new = [[] for _ in range(self.W)]
            for c in range(self.W):
                bits = list(cols[c])
                while len(bits) >= 3:
                    s, co = self.fa(bits.pop(), bits.pop(), bits.pop())
                    new[c].append(s)
                    if c + 1 < self.W: new[c + 1].append(co)
                if len(bits) == 2:
                    s, co = self.ha(bits.pop(), bits.pop())
                    new[c].append(s)
                    if c + 1 < self.W: new[c + 1].append(co)
                new[c].extend(bits)
            cols = new
            if lvl in self.pipe:
                cols = [[self.reg(v) for v in col] for col in cols]
            lvl += 1
        S = {c: col[0] for c, col in enumerate(cols) if len(col) > 0}
        C = {c: col[1] for c, col in enumerate(cols) if len(col) > 1}
        return S, C

    def dadda(self, a, b):
        cols = self.booth_cols(a, b) if self.arch == "booth" else self.bw_cols(a, b)
        seq = [2]
        while seq[-1] < max(len(c) for c in cols):
            seq.append(seq[-1] * 3 // 2)
        targets = [d for d in reversed(seq) if d < max(len(c) for c in cols)]
        for lvl, d in enumerate(targets):
            new = [[] for _ in range(self.W)]
            for c in range(self.W):
                bits = list(cols[c])
                while len(bits) + len(new[c]) > d:
                    excess = len(bits) + len(new[c]) - d
                    if excess >= 2 and len(bits) >= 3:
                        s, co = self.fa(bits.pop(), bits.pop(), bits.pop())
                    else:
                        s, co = self.ha(bits.pop(), bits.pop())
                    new[c].append(s)
                    if c + 1 < self.W:
                        new[c + 1].append(co)
                new[c].extend(bits)
            cols = new
            if lvl in self.pipe:
                cols = [[self.reg(v) for v in col] for col in cols]
        S = {c: col[0] for c, col in enumerate(cols) if len(col) > 0}
        C = {c: col[1] for c, col in enumerate(cols) if len(col) > 1}
        return S, C

    def ks_adder(self, S, C):
        """Kogge-Stone parallel-prefix vector-merging adder (log2(W) prefix levels)."""
        W = self.W
        g, p = [], []
        for c in range(W):
            x, y = S.get(c), C.get(c)
            if x and y:
                g.append(self.w(f"{x} & {y}", 1, (x, y))); p.append(self.w(f"{x} ^ {y}", 1, (x, y)))
            else:
                g.append("1'b0"); p.append(x or y or "1'b0")
        p0 = list(p)
        G, P = list(g), list(p)
        d = 1
        while d < W:
            nG, nP = list(G), list(P)
            for i in range(d, W):
                if G[i - d] != "1'b0":
                    nG[i] = self.w(f"{G[i]} | ({P[i]} & {G[i - d]})", 2, (G[i], P[i], G[i - d]))
                nP[i] = self.w(f"{P[i]} & {P[i - d]}", 1, (P[i], P[i - d]))
            G, P = nG, nP
            d *= 2
        out = [p0[0]] + [self.w(f"{p0[i]} ^ {G[i - 1]}", 1, (p0[i], G[i - 1])) if G[i - 1] != "1'b0" else p0[i] for i in range(1, W)]
        return out

    def final_adder(self, S, C):
        """Ripple-carry vector-merging adder; optional register split at the midpoint."""
        if self.cpa == "ks":
            return self.ks_adder(S, C)
        out, carry = [], None
        split = self.W // 2 if "F" in self.pipe else None
        for c in range(self.W):
            if c == split:
                out = [self.reg(x) for x in out]
                S = {k: self.reg(v) for k, v in S.items() if k >= c}
                C = {k: self.reg(v) for k, v in C.items() if k >= c}
                carry = self.reg(carry) if carry else None
            bits = [x for x in (S.get(c), C.get(c), carry) if x]
            if len(bits) == 3:
                s, carry = self.fa(*bits)
            elif len(bits) == 2:
                s, carry = self.ha(*bits)
            elif len(bits) == 1:
                s, carry = bits[0], None
            else:
                s, carry = "1'b0", None
            out.append(s)
        return out

    def emit(self):
        n = self.n
        a = [f"a[{j}]" for j in range(n)]
        b = [f"b[{i}]" for i in range(n)]
        red = {"array": self.array, "dadda": self.dadda, "booth": self.dadda, "wallace": self.wallace}
        S, C = red[self.arch](a, b)
        self.tree_depth = max(self.depth.get(x, 0) for x in list(S.values()) + list(C.values()))
        p = self.final_adder(S, C)
        self.out_bits = p
        seq = bool(self.regs)
        ports = "input wire clk, " if seq else ""
        hdr = [f"// Auto-generated by gen_mult.py: {self.arch}{' (radix-4 Booth + Dadda)' if self.arch == 'booth' else ''} + {self.cpa.upper()} final adder, {('signed' if self.arch == 'booth' else 'signed Baugh-Wooley') if self.signed else 'unsigned'}, "
               f"{n}x{n}, pipeline={sorted(map(str, self.pipe))}",
               f"module {self.name} ({ports}input wire [{n-1}:0] a, input wire [{n-1}:0] b, output wire [{self.W-1}:0] p);"]
        body = list(self.lines)
        if seq:
            body += [f"  reg {r};" for r, _ in self.regs]
            body.append("  always @(posedge clk) begin")
            body += [f"    {r} <= {s};" for r, s in self.regs]
            body.append("  end")
        body.append(f"  assign p = {{{', '.join(reversed(p))}}};")
        return "\n".join(hdr + body + ["endmodule", ""])

    def stats(self):
        """Analytical structure: adder counts, partial-product gates, max stage logic depth."""
        depth_out = max(self.depth.get(x, 0) for x in self.out_bits)
        return dict(fa=self.nfa, ha=self.nha, regs=len(self.regs), tree_depth=self.tree_depth,
                    depth=max([depth_out] + self.stage_depths))

    def latency(self):
        # number of distinct register stages = number of pipeline cut points used
        return len(self.pipe)

if __name__ == "__main__":
    name, n, arch, sg = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4] == "s"
    pipe = []
    if len(sys.argv) > 5:
        pipe = [x if x == "F" else int(x) for x in sys.argv[5].split(",")]
    g = Gen(name, n, arch, sg, pipe)
    sys.stdout.write(g.emit())
