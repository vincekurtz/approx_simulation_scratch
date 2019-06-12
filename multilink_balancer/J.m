function J_sym = J(in1)
%J
%    J_SYM = J(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Jun-2019 10:20:13

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
t2 = theta1+theta2+theta3;
t3 = sin(t2);
t4 = theta1+theta2+theta4;
t5 = sin(t4);
t6 = theta1+theta2;
t7 = sin(t6);
t8 = cos(t2);
t9 = t8./8.0;
t10 = cos(t4);
t11 = t10./8.0;
t12 = cos(t6);
t13 = t12.*(5.0./8.0);
J_sym = reshape([t3.*(-1.0./8.0)-t5./8.0-t7.*(5.0./8.0)-sin(theta1).*(7.0./8.0),t9+t11+t13+cos(theta1).*(7.0./8.0),t3.*(-1.0./8.0)-t5./8.0-t7.*(5.0./8.0),t9+t11+t13,t3.*(-1.0./8.0),t9,t5.*(-1.0./8.0),t11],[2,4]);
