function J_sym = J(in1)
%J
%    J_SYM = J(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    21-May-2019 16:39:01

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
t2 = theta1+theta2+theta3;
t3 = sin(t2);
t4 = theta1+theta2;
t5 = sin(t4);
t6 = cos(t2);
t7 = t6./6.0;
t8 = cos(t4);
t9 = t8./2.0;
J_sym = reshape([t3.*(-1.0./6.0)-t5./2.0-sin(theta1).*(5.0./6.0),t7+t9+cos(theta1).*(5.0./6.0),t3.*(-1.0./6.0)-t5./2.0,t7+t9,t3.*(-1.0./6.0),t7],[2,3]);