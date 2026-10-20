// Behavioural reference: synthesis-tool default implementation of a signed 4x4 multiply.
module behav_mult4 (input wire [3:0] a, input wire [3:0] b, output wire [7:0] p);
  assign p = $signed(a) * $signed(b);
endmodule
