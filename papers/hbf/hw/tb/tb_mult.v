// Exhaustive testbench for generated multipliers.
// Defines: DUT (module), N (width), SIGNED (1/0), LAT (pipeline latency, 0 = combinational)
`timescale 1ns/1ps
module tb;
  localparam N = `N, W = 2*N, LAT = `LAT, EXH = (2*N <= 16), TOTAL = EXH ? (1 << (2*N)) : 300000;
  integer seed = 1;
  reg clk = 0;
  reg [N-1:0] a, b;
  wire [W-1:0] p;
  reg [W-1:0] expq [0:LAT];
  integer k, i, errors = 0, checked = 0;
`ifdef PIPE
  `DUT dut (.clk(clk), .a(a), .b(b), .p(p));
`else
  `DUT dut (.a(a), .b(b), .p(p));
`endif
  function [W-1:0] mref(input [N-1:0] x, input [N-1:0] y);
    if (`SIGNED) mref = $signed(x) * $signed(y); else mref = x * y;
  endfunction
  initial begin
    for (k = 0; k < TOTAL + LAT; k = k + 1) begin
      if (EXH) begin a = k[N-1:0]; b = k[2*N-1:N]; end
      else if (k < 16) begin   // corners: 0, 1, -1, min, max combinations
        a = (k[1:0] == 0) ? 0 : (k[1:0] == 1) ? {1'b0, {(N-1){1'b1}}} : (k[1:0] == 2) ? {1'b1, {(N-1){1'b0}}} : {N{1'b1}};
        b = (k[3:2] == 0) ? 1 : (k[3:2] == 1) ? {1'b0, {(N-1){1'b1}}} : (k[3:2] == 2) ? {1'b1, {(N-1){1'b0}}} : {N{1'b1}};
      end else begin a = $random(seed); b = $random(seed); end
      if (LAT == 0) begin
        #1;
        if (p !== mref(a, b)) begin errors = errors + 1;
          if (errors < 5) $display("MISMATCH a=%0d b=%0d p=%h exp=%h", a, b, p, mref(a, b)); end
        checked = checked + 1;
      end else begin
        for (i = LAT; i > 0; i = i - 1) expq[i] = expq[i-1];
        expq[0] = mref(a, b);
        #1 clk = 1; #1 clk = 0;
        if (k >= LAT-1 && k < TOTAL+LAT-1) begin
          if (p !== expq[LAT-1]) begin errors = errors + 1;
            if (errors < 5) $display("MISMATCH k=%0d p=%h exp=%h", k, p, expq[LAT-1]); end
          checked = checked + 1;
        end
      end
    end
    $display("RESULT: %0d vectors checked, %0d errors -> %s", checked, errors, errors ? "FAIL" : "PASS");
    $finish;
  end
endmodule
