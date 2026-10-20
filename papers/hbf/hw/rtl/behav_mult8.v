// Behavioural reference: synthesis-tool default implementation of a signed 8x8 multiply.
module behav_mult8 (input wire [7:0] a, input wire [7:0] b, output wire [15:0] p);
  assign p = $signed(a) * $signed(b);
endmodule
