function vpc_sym = pd_com(in1,in2)
%PD_COM
%    VPC_SYM = PD_COM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    18-Jun-2019 14:53:03

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
theta_dot1 = in2(1,:);
theta_dot2 = in2(2,:);
theta_dot3 = in2(3,:);
theta_dot4 = in2(4,:);
t2 = theta1+theta2+theta3;
t3 = sin(t2);
t4 = theta1+theta2+theta3+theta4;
t5 = sin(t4);
t6 = theta1+theta2;
t7 = sin(t6);
t8 = cos(t4);
t9 = cos(t6);
t10 = cos(t2);
vpc_sym = [t3.*theta_dot1.*(-3.0./5.0)-t3.*theta_dot2.*(3.0./5.0)-t3.*theta_dot3.*(3.0./5.0)-(t5.*theta_dot1)./1.0e1-(t5.*theta_dot2)./1.0e1-(t5.*theta_dot3)./1.0e1-t7.*theta_dot1.*(7.0./1.0e1)-(t5.*theta_dot4)./1.0e1-t7.*theta_dot2.*(7.0./1.0e1)-theta_dot1.*sin(theta1).*(9.0./1.0e1);(t8.*theta_dot1)./1.0e1+(t8.*theta_dot2)./1.0e1+t9.*theta_dot1.*(7.0./1.0e1)+(t8.*theta_dot3)./1.0e1+t9.*theta_dot2.*(7.0./1.0e1)+t10.*theta_dot1.*(3.0./5.0)+(t8.*theta_dot4)./1.0e1+t10.*theta_dot2.*(3.0./5.0)+t10.*theta_dot3.*(3.0./5.0)+theta_dot1.*cos(theta1).*(9.0./1.0e1)];
