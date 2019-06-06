function pc_sym = p_com(in1)
%P_COM
%    PC_SYM = P_COM(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    06-Jun-2019 10:27:28

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
t2 = theta1+theta2+theta3;
t3 = theta1+theta2+theta4;
t4 = theta1+theta2;
pc_sym = [cos(t2)./8.0+cos(t3)./8.0+cos(t4).*(5.0./8.0)+cos(theta1).*(7.0./8.0);sin(t2)./8.0+sin(t3)./8.0+sin(t4).*(5.0./8.0)+sin(theta1).*(7.0./8.0)];
