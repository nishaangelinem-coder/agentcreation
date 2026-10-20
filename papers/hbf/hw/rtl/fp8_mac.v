// FP8 (E4M3) multiply-accumulate with exact fixed-point accumulation.
// Product of two E4M3 values = (m1*m2) * 2^(E1+E2-20), E = max(exp,1), m = {exp!=0, mant}.
// Expressed in units of 2^-18 the product is (m1*m2) << (E1+E2-2), shift 0..28, so a
// 37-bit signed term is exact; a 48-bit accumulator gives 2^11 headroom (>= 2048 terms
// at full scale). The E4M3 NaN code (S.1111.111) is treated as the finite value 480.
// Pipeline: S1 decode + 4x4 unsigned array multiply + exponent add | S2 align | S3 accumulate.
module fp8_mac #(parameter ACC_W = 48) (
  input  wire              clk,
  input  wire              rst,
  input  wire              clr,      // start a new dot product with this input
  input  wire              valid,
  input  wire [7:0]        w,        // weight (streamed from HBF page buffer)
  input  wire [7:0]        x,        // activation (broadcast)
  output reg  [ACC_W-1:0]  acc,
  output reg  signed [36:0] term     // aligned signed product (exposed for verification)
);
  wire [3:0] mw = {|w[6:3], w[2:0]};
  wire [3:0] mx = {|x[6:3], x[2:0]};
  wire [3:0] ew = (w[6:3] == 4'd0) ? 4'd1 : w[6:3];
  wire [3:0] ex = (x[6:3] == 4'd0) ? 4'd1 : x[6:3];
  wire [7:0] mp;
  u_array4 mul (.a(mw), .b(mx), .p(mp));

  // stage 1
  reg [7:0] s1_mp; reg [4:0] s1_sh; reg s1_sign, s1_v, s1_clr;
  always @(posedge clk) begin
    s1_mp <= mp; s1_sh <= ew + ex - 5'd2; s1_sign <= w[7] ^ x[7];
    s1_v <= valid & ~rst; s1_clr <= clr;
  end
  // stage 2: barrel shift and conditional negate
  wire [35:0] mag = {28'd0, s1_mp} << s1_sh;
  reg s2_v, s2_clr;
  always @(posedge clk) begin
    term <= s1_sign ? -$signed({1'b0, mag}) : $signed({1'b0, mag});
    s2_v <= s1_v & ~rst; s2_clr <= s1_clr;
  end
  // stage 3: accumulate
  wire [ACC_W-1:0] term_x = {{(ACC_W-37){term[36]}}, term};
  always @(posedge clk) begin
    if (rst) acc <= 0;
    else if (s2_v) acc <= s2_clr ? term_x : acc + term_x;
  end
endmodule
