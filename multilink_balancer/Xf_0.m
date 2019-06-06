function com_Xf_o = Xf_0(in1)
%XF_0
%    COM_XF_O = XF_0(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    06-Jun-2019 10:29:20

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
t2 = theta1+theta2+theta3;
t3 = theta1+theta2+theta4;
t4 = theta1+theta2;
t5 = sin(t2);
t6 = sin(t3);
t7 = sin(t4);
t8 = sin(theta1);
t9 = cos(t2);
t10 = t9./8.0;
t11 = cos(t3);
t12 = t11./8.0;
t13 = cos(t4);
t14 = t13.*(5.0./8.0);
t15 = cos(theta1);
t16 = t15.*(7.0./8.0);
com_Xf_o = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,t5./8.0+t6./8.0+t7.*(5.0./8.0)+t8.*(7.0./8.0),1.0,0.0,0.0,0.0,0.0,-t10-t12-t14-t16,0.0,1.0,0.0,t5.*(-1.0./8.0)-t6./8.0-t7.*(5.0./8.0)-t8.*(7.0./8.0),t10+t12+t14+t16,0.0,0.0,0.0,1.0],[6,6]);
