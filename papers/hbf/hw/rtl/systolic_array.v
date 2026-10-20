// ROWS x COLS output-stationary systolic array for the HBF base die.
// Each column is fed by one NAND channel lane (one weight byte per cycle);
// each row holds one sequence of the decode batch (activations from an on-die SRAM).
// PE(r,c) accumulates sum_k A[r][k] * W[k][c]. Inputs must be skewed by the caller:
// row r is delayed r cycles and column c is delayed c cycles.
module systolic_array #(parameter ROWS = 4, COLS = 8, ACC_W = 32) (
  input  wire                       clk,
  input  wire                       clr,
  input  wire [ROWS*8-1:0]          a_west,
  input  wire [COLS*8-1:0]          w_north,
  output wire [ROWS*COLS*ACC_W-1:0] acc_flat
);
  wire [7:0] ah [0:ROWS-1][0:COLS];
  wire [7:0] wv [0:ROWS][0:COLS-1];
  genvar r, c;
  generate
    for (r = 0; r < ROWS; r = r + 1) begin : g_west
      assign ah[r][0] = a_west[r*8 +: 8];
    end
    for (c = 0; c < COLS; c = c + 1) begin : g_north
      assign wv[0][c] = w_north[c*8 +: 8];
    end
    for (r = 0; r < ROWS; r = r + 1) begin : g_r
      for (c = 0; c < COLS; c = c + 1) begin : g_c
        wire signed [ACC_W-1:0] acc;
        int8_pe #(.ACC_W(ACC_W)) pe (.clk(clk), .clr(clr), .a_in(ah[r][c]), .w_in(wv[r][c]),
                                     .a_out(ah[r][c+1]), .w_out(wv[r+1][c]), .acc(acc));
        assign acc_flat[(r*COLS + c)*ACC_W +: ACC_W] = acc;
      end
    end
  endgenerate
endmodule
