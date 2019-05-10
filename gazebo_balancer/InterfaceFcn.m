function tau_sym = InterfaceFcn(u_lip,in2,in3)
%INTERFACEFCN
%    TAU_SYM = INTERFACEFCN(U_LIP,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2019 16:28:56

theta1 = in3(1,:);
theta2 = in3(2,:);
theta1_dot = in3(3,:);
theta2_dot = in3(4,:);
x_lip1 = in2(1,:);
x_lip2 = in2(2,:);
x_lip3 = in2(3,:);
x_lip4 = in2(4,:);
t2 = cos(theta1);
t3 = cos(theta2);
t4 = sin(theta2);
t5 = theta1+theta2;
t6 = sin(theta1);
t7 = sin(t5);
t8 = cos(t5);
t9 = t2.*t7;
t14 = t6.*t8;
t10 = t9-t14;
t11 = 1.0./t10;
t12 = t3./2.0;
t13 = theta1_dot+theta2_dot;
t15 = t12+3.345833333333333e-1;
t16 = t4.^2;
t17 = t3+1.0./2.0;
t18 = t3.*t17;
t19 = t12+t16+t18+6.691666666666667e-1;
t20 = t8./4.0;
t21 = t2.*(3.0./4.0);
t22 = sqrt(3.0);
t23 = t7./4.0;
t24 = t6.*(3.0./4.0);
t25 = t2.*t3.*(9.81e2./2.0e2);
t26 = t2.*3.0;
t27 = t8+t26;
t28 = x_lip1.*(3.52e2./2.5e1);
t29 = (t8.*t13)./4.0;
t30 = t2.*theta1_dot.*(3.0./4.0);
t31 = t29+t30;
t32 = t31.*theta1_dot;
t33 = (t7.*theta2_dot)./4.0;
t34 = t23+t24;
t35 = t34.*theta1_dot;
t36 = t33+t35+x_lip2;
t37 = t22.*t36;
t38 = (t8.*t13.*theta2_dot)./4.0;
t39 = -t20-t21+t28+t32+t37+t38-u_lip.*(3.27e2./2.5e1);
t40 = t6.*3.0;
t41 = t7+t40;
t42 = (t7.*t13)./4.0;
t43 = t6.*theta1_dot.*(3.0./4.0);
t44 = t42+t43;
t45 = t20+t21;
t46 = t45.*theta1_dot;
t47 = (t8.*theta2_dot)./4.0;
t48 = t46+t47-x_lip4;
t49 = t22.*t48;
t50 = (t4.*t13.*theta1_dot)./2.0;
tau_sym = [t2.*(9.81e2./2.0e2)+t25+t50+t4.*(t2.*t4.*(9.81e2./1.0e2)+t3.*t6.*(9.81e2./1.0e2)-t13.*(theta1_dot./2.0+theta2_dot./2.0)-t3.*t13.*theta1_dot+t3.*theta1_dot.*theta2_dot)+t3.*(t2.*t3.*(9.81e2./1.0e2)-t4.*t6.*(9.81e2./1.0e2)+t4.*t13.*theta1_dot-t4.*theta1_dot.*theta2_dot)-t4.*t6.*(9.81e2./2.0e2)+t39.*(t8.*t11.*t19.*(4.0./3.0)-t11.*t15.*t27.*(4.0./3.0))-(t7.*t11.*t19.*(4.0./3.0)-t11.*t15.*t41.*(4.0./3.0)).*(t23+t24+t49-x_lip3-t44.*theta1_dot-(t7.*t13.*theta2_dot)./4.0)-(t4.*theta1_dot.*theta2_dot)./2.0;t25+t50-t39.*(t11.*t27.*4.461111111111111e-1-t8.*t11.*t15.*(4.0./3.0))-t4.*t6.*(9.81e2./2.0e2)+(t11.*t41.*4.461111111111111e-1-t7.*t11.*t15.*(4.0./3.0)).*(t23+t24+t49-x_lip3-t44.*theta1_dot-(t7.*t13.*theta2_dot)./4.0)-(t4.*theta1_dot.*theta2_dot)./2.0];
