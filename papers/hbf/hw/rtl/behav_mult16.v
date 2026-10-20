// Behavioural reference: synthesis-tool default implementation of a signed 16x16 multiply.
module behav_mult16 (input wire [15:0] a, input wire [15:0] b, output wire [31:0] p);
  assign p = $signed(a) * $signed(b);
endmodule
