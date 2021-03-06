function f_sym = BalancerDynamics(in1,in2)
%BALANCERDYNAMICS
%    F_SYM = BALANCERDYNAMICS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    26-Apr-2019 16:53:25

tau1 = in2(1,:);
tau2 = in2(2,:);
theta1 = in1(1,:);
theta2 = in1(2,:);
theta1_dot = in1(3,:);
theta2_dot = in1(4,:);
t2 = cos(theta2);
t3 = sin(theta2);
t4 = theta1_dot+theta2_dot;
t5 = t2.^2;
t6 = t5.*(2.03e2./8.03e2);
t7 = t3.^2;
t8 = cos(theta1);
t9 = t6+t7+3.345833333333333e-1;
t10 = 1.0./t9;
t11 = theta1_dot./2.0;
t12 = theta2_dot./2.0;
t13 = t11+t12;
t14 = t4.*t13;
t15 = t2.*t4.*theta1_dot;
t26 = t2.*theta1_dot.*theta2_dot;
t16 = t14+t15-t26;
t17 = t6+t7+1.0./2.0;
t18 = t8.*t17.*(9.81e2./1.0e2);
t19 = tau2.*1.49439601494396;
t20 = t3.*t4.*theta1_dot.*(2.03e2./8.03e2);
t28 = t3.*theta1_dot.*theta2_dot.*(2.03e2./8.03e2);
t21 = t19+t20-t28;
t22 = t2.*t21;
t23 = sin(theta1);
t24 = t2.*t3.*t23.*7.330012453300125;
t27 = t3.*t16;
t25 = t18+t22+t24-t27-tau1+tau2;
f_sym = [theta1_dot;theta2_dot;-t10.*t25;tau2.*2.98879202988792-t2.*t8.*1.466002490660025e1+t3.*t23.*1.466002490660025e1+t10.*t25+t2.*t10.*t25.*1.49439601494396-t3.*t4.*theta1_dot.*1.49439601494396+t3.*theta1_dot.*theta2_dot.*1.49439601494396];
