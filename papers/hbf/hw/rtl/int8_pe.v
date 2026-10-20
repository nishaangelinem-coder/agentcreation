// Output-stationary INT8 processing element: activations move right, weights move down.
// Product is registered once (multiplier | accumulate) to keep the cycle to one array multiply.
module int8_pe #(parameter ACC_W = 32) (
  input  wire                    clk,
  input  wire                    clr,
  input  wire signed [7:0]       a_in,
  input  wire signed [7:0]       w_in,
  output reg  signed [7:0]       a_out,
  output reg  signed [7:0]       w_out,
  output reg  signed [ACC_W-1:0] acc
);
  wire [15:0] p;
  bw_array8 mul (.a(a_in), .b(w_in), .p(p));
  reg signed [15:0] p_r; reg clr_r;
  always @(posedge clk) begin
    a_out <= a_in; w_out <= w_in;
    p_r <= p; clr_r <= clr;
    acc <= clr_r ? {ACC_W{1'b0}} : acc + {{(ACC_W-16){p_r[15]}}, p_r};
  end
endmodule
