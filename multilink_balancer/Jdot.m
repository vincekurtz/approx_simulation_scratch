function Jdot_sym = Jdot(in1,in2)
%JDOT
%    JDOT_SYM = JDOT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    21-May-2019 16:39:01

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta_dot1 = in2(1,:);
theta_dot2 = in2(2,:);
theta_dot3 = in2(3,:);
t2 = theta1+theta2;
t3 = cos(t2);
t4 = theta1+theta2+theta3;
t5 = cos(t4);
t6 = sin(t4);
t7 = sin(t2);
Jdot_sym = reshape([t3.*theta_dot1.*(-1.0./2.0)-(t3.*theta_dot2)./2.0-(t5.*theta_dot1)./6.0-(t5.*theta_dot2)./6.0-(t5.*theta_dot3)./6.0-theta_dot1.*cos(theta1).*(5.0./6.0),t6.*theta_dot1.*(-1.0./6.0)-(t6.*theta_dot2)./6.0-(t7.*theta_dot1)./2.0-(t6.*theta_dot3)./6.0-(t7.*theta_dot2)./2.0-theta_dot1.*sin(theta1).*(5.0./6.0),t3.*theta_dot1.*(-1.0./2.0)-(t3.*theta_dot2)./2.0-(t5.*theta_dot1)./6.0-(t5.*theta_dot2)./6.0-(t5.*theta_dot3)./6.0,t6.*theta_dot1.*(-1.0./6.0)-(t6.*theta_dot2)./6.0-(t7.*theta_dot1)./2.0-(t6.*theta_dot3)./6.0-(t7.*theta_dot2)./2.0,t5.*theta_dot1.*(-1.0./6.0)-(t5.*theta_dot2)./6.0-(t5.*theta_dot3)./6.0,t6.*theta_dot1.*(-1.0./6.0)-(t6.*theta_dot2)./6.0-(t6.*theta_dot3)./6.0],[2,3]);
