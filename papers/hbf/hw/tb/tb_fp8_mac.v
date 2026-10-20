`timescale 1ns/1ps
module tb;
  reg clk = 0, rst = 1, clr = 0, valid = 0;
  reg [7:0] w = 0, x = 0;
  wire [47:0] acc; wire signed [36:0] term;
  fp8_mac dut (.clk(clk), .rst(rst), .clr(clr), .valid(valid), .w(w), .x(x), .acc(acc), .term(term));
  reg [36:0] prod_mem [0:65535];
  reg [47:0] dot_mem [0:300000];
  integer k, errors = 0, idx, L, t, n, dots = 0;
  task tick; begin #1 clk = 1; #1 clk = 0; end endtask
  initial begin
    $readmemh("build/fp8_products.hex", prod_mem);
    $readmemh("build/fp8_dots.hex", dot_mem);
    tick; rst = 0;
    // exhaustive product check (term is valid 2 cycles after inputs)
    for (k = 0; k < 65536 + 2; k = k + 1) begin
      w = k[7:0]; x = k[15:8]; valid = 1; tick;
      if (k >= 1 && k <= 65536 && term !== prod_mem[k-1]) begin
        errors = errors + 1; if (errors < 5) $display("TERM MISMATCH k=%0d got=%h exp=%h", k-1, term, prod_mem[k-1]); end
    end
    // dot products
    idx = 0;
    for (t = 0; t < 200; t = t + 1) begin
      L = dot_mem[idx]; idx = idx + 1;
      for (n = 0; n < L; n = n + 1) begin
        {w, x} = dot_mem[idx][15:0]; idx = idx + 1; valid = 1; clr = (n == 0); tick;
      end
      valid = 0; clr = 0; tick; tick; tick;
      if (acc !== dot_mem[idx]) begin errors = errors + 1;
        if (errors < 5) $display("DOT MISMATCH t=%0d L=%0d got=%h exp=%h", t, L, acc, dot_mem[idx]); end
      idx = idx + 1; dots = dots + 1;
    end
    $display("RESULT: 65536 products + %0d dot products checked, %0d errors -> %s", dots, errors, errors ? "FAIL" : "PASS");
    $finish;
  end
endmodule
