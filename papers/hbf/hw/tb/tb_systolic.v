`timescale 1ns/1ps
module tb;
  localparam R = 4, C = 8, K = 64, AW = 32, TRIALS = 50;
  reg clk = 0, clr = 0;
  reg [R*8-1:0] a_west; reg [C*8-1:0] w_north;
  wire [R*C*AW-1:0] acc_flat;
  systolic_array #(.ROWS(R), .COLS(C), .ACC_W(AW)) dut (.clk(clk), .clr(clr), .a_west(a_west), .w_north(w_north), .acc_flat(acc_flat));
  reg signed [7:0] A [0:R-1][0:K-1];
  reg signed [7:0] Wt [0:K-1][0:C-1];
  integer r, c, k, t, cyc, errors = 0, checks = 0, seed = 11;
  reg signed [31:0] exp_v;
  task tick; begin #1 clk = 1; #1 clk = 0; end endtask
  initial begin
    for (t = 0; t < TRIALS; t = t + 1) begin
      for (r = 0; r < R; r = r + 1) for (k = 0; k < K; k = k + 1) A[r][k] = (t == 0) ? -128 : $random(seed);
      for (k = 0; k < K; k = k + 1) for (c = 0; c < C; c = c + 1) Wt[k][c] = (t == 0) ? -128 : $random(seed);
      a_west = 0; w_north = 0; clr = 1; for (cyc = 0; cyc < R + C + 2; cyc = cyc + 1) tick; clr = 0;
      for (cyc = 0; cyc < K + R + C + 4; cyc = cyc + 1) begin
        for (r = 0; r < R; r = r + 1) a_west[r*8 +: 8] = (cyc - r >= 0 && cyc - r < K) ? A[r][cyc - r] : 8'd0;
        for (c = 0; c < C; c = c + 1) w_north[c*8 +: 8] = (cyc - c >= 0 && cyc - c < K) ? Wt[cyc - c][c] : 8'd0;
        tick;
      end
      for (r = 0; r < R; r = r + 1) for (c = 0; c < C; c = c + 1) begin
        exp_v = 0; for (k = 0; k < K; k = k + 1) exp_v = exp_v + A[r][k] * Wt[k][c];
        if ($signed(acc_flat[(r*C+c)*AW +: AW]) !== exp_v) begin errors = errors + 1;
          if (errors < 5) $display("MISMATCH t=%0d r=%0d c=%0d got=%0d exp=%0d", t, r, c, $signed(acc_flat[(r*C+c)*AW +: AW]), exp_v); end
        checks = checks + 1;
      end
    end
    $display("RESULT: %0d trials, %0d outputs checked, %0d errors -> %s", TRIALS, checks, errors, errors ? "FAIL" : "PASS");
    $finish;
  end
endmodule
