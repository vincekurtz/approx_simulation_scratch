function h_o = h_o(in1,in2)
%H_O
%    H_O = H_O(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Jun-2019 10:22:09

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
theta_dot1 = in2(1,:);
theta_dot2 = in2(2,:);
theta_dot3 = in2(3,:);
theta_dot4 = in2(4,:);
t2 = sin(theta3);
t3 = cos(theta2);
t4 = sin(theta4);
t5 = theta_dot1+theta_dot2;
t6 = sin(theta2);
t7 = cos(theta3);
t8 = theta_dot1./2.0;
t9 = theta_dot2./2.0;
t10 = cos(theta4);
t11 = t2.*t5;
t12 = t2.*t3.*theta_dot1;
t13 = t6.*t7.*theta_dot1;
t14 = t11+t12+t13;
t15 = t2.*t14;
t16 = t4.*t5;
t17 = t3.*t4.*theta_dot1;
t18 = t6.*t10.*theta_dot1;
t19 = t16+t17+t18;
t20 = t4.*t19;
t21 = theta_dot3./2.0;
t22 = t5.*t7;
t23 = t3.*t7.*theta_dot1;
t29 = t2.*t6.*theta_dot1;
t24 = t8+t9+t21+t22+t23-t29;
t25 = theta_dot4./2.0;
t26 = t5.*t10;
t27 = t3.*t10.*theta_dot1;
t31 = t4.*t6.*theta_dot1;
t28 = t8+t9+t25+t26+t27-t31;
t30 = t7.*t24;
t32 = t10.*t28;
t33 = t3.*theta_dot1;
t34 = t8+t9+t15+t20+t30+t32+t33;
t35 = t3.*t34;
t36 = t7.*t14;
t37 = t10.*t19;
t38 = t6.*theta_dot1;
t41 = t2.*t24;
t42 = t4.*t28;
t39 = t36+t37+t38-t41-t42;
t40 = t6.*t39;
t43 = cos(theta1);
t44 = t8+t35+t40;
t45 = sin(theta1);
t46 = t6.*t34;
t47 = t46-t3.*t39;
h_o = [0.0;0.0;t15+t20+t30+t32+t35+t40+theta_dot1.*(8.03e2./6.0e2)+theta_dot2.*(8.03e2./8.0e2)+theta_dot3.*3.345833333333333e-1+theta_dot4.*3.345833333333333e-1+(t5.*t7)./2.0+(t5.*t10)./2.0+(t3.*theta_dot1)./2.0-(t2.*t6.*theta_dot1)./2.0+(t3.*t7.*theta_dot1)./2.0-(t4.*t6.*theta_dot1)./2.0+(t3.*t10.*theta_dot1)./2.0;-t44.*t45-t43.*t47;t43.*t44-t45.*t47;0.0];
