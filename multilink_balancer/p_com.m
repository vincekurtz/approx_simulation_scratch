function pc_sym = p_com(in1)
%P_COM
%    PC_SYM = P_COM(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    22-May-2019 09:14:05

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
t2 = theta1+theta2+theta3;
t3 = theta1+theta2;
pc_sym = [cos(t2)./6.0+cos(t3)./2.0+cos(theta1).*(5.0./6.0);sin(t2)./6.0+sin(t3)./2.0+sin(theta1).*(5.0./6.0)];
