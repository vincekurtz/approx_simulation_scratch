function A_cwc_sym = A_cwc(in1)
%A_CWC
%    A_CWC_SYM = A_CWC(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Jun-2019 10:22:16

p_com1 = in1(1,:);
p_com2 = in1(2,:);
A_cwc_sym = reshape([0.0,0.0,0.0,1.0,-1.0,0.0,1.0,-1.0,-p_com2,p_com2,-1.0,-1.0./5.0,-1.0./5.0,p_com1-1.0./2.0,-p_com1-1.0./2.0],[5,3]);
