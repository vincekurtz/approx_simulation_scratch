function C_sym = C(in1,in2)
%C
%    C_SYM = C(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    24-May-2019 10:27:12

theta1 = in1(1,:);
theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
theta_dot1 = in2(1,:);
theta_dot2 = in2(2,:);
theta_dot3 = in2(3,:);
theta_dot4 = in2(4,:);
t2 = cos(theta1);
t3 = sin(theta2);
t4 = cos(theta2);
t5 = sin(theta1);
t6 = t3.*t5.*(9.81e2./1.0e2);
t7 = t3.*theta_dot1.*theta_dot2;
t18 = t2.*t4.*(9.81e2./1.0e2);
t8 = t6+t7-t18;
t9 = t2.*t3.*(9.81e2./1.0e2);
t10 = t4.*t5.*(9.81e2./1.0e2);
t11 = t4.*theta_dot1.*theta_dot2;
t12 = t9+t10+t11;
t13 = sin(theta3);
t14 = cos(theta3);
t15 = sin(theta4);
t16 = theta_dot1+theta_dot2;
t17 = cos(theta4);
t19 = theta_dot1+theta_dot2+theta_dot3;
t20 = t14.*t16;
t21 = t4.*t14.*theta_dot1;
t38 = t3.*t13.*theta_dot1;
t22 = t20+t21-t38;
t23 = theta_dot1+theta_dot2+theta_dot4;
t24 = t16.*t17;
t25 = t4.*t17.*theta_dot1;
t45 = t3.*t15.*theta_dot1;
t26 = t24+t25-t45;
t27 = theta_dot1./2.0;
t28 = theta_dot2./2.0;
t29 = t13.*t16;
t30 = t4.*t13.*theta_dot1;
t31 = t3.*t14.*theta_dot1;
t32 = t29+t30+t31;
t33 = t15.*t16;
t34 = t4.*t15.*theta_dot1;
t35 = t3.*t17.*theta_dot1;
t36 = t33+t34+t35;
t37 = t8.*t13;
t39 = t19.*t22;
t40 = theta_dot3./2.0;
t41 = t27+t28+t40;
t42 = t19.*t41;
t59 = t22.*theta_dot3;
t60 = t12.*t14;
t43 = t37+t39+t42-t59-t60;
t44 = t8.*t15;
t46 = t23.*t26;
t47 = theta_dot4./2.0;
t48 = t27+t28+t47;
t49 = t23.*t48;
t61 = t26.*theta_dot4;
t62 = t12.*t17;
t50 = t44+t46+t49-t61-t62;
t51 = t8.*t14;
t52 = t12.*t13;
t53 = t32.*theta_dot3;
t63 = t19.*t32;
t54 = t51+t52+t53-t63;
t55 = t8.*t17;
t56 = t12.*t15;
t57 = t36.*theta_dot4;
t64 = t23.*t36;
t58 = t55+t56+t57-t64;
t65 = t2.*t4.*(9.81e2./2.0e2);
t66 = (t13.*t16)./2.0;
t67 = (t4.*t13.*theta_dot1)./2.0;
t68 = (t3.*t14.*theta_dot1)./2.0;
t69 = t66+t67+t68;
t70 = t19.*t69;
t71 = (t15.*t16)./2.0;
t72 = (t4.*t15.*theta_dot1)./2.0;
t73 = (t3.*t17.*theta_dot1)./2.0;
t74 = t71+t72+t73;
t75 = t23.*t74;
t76 = t13.*t43;
t77 = t15.*t50;
t78 = t14.*t54;
t79 = t17.*t58;
t80 = (t3.*t16.*theta_dot1)./2.0;
C_sym = [t2.*(9.81e2./2.0e2)+t65+t70+t75+t80-t3.*t5.*(9.81e2./2.0e2)-(t8.*t14)./2.0-(t8.*t17)./2.0-(t12.*t13)./2.0-(t12.*t15)./2.0-t13.*t43-t15.*t50-t14.*t54-t17.*t58-(t32.*theta_dot3)./2.0-(t36.*theta_dot4)./2.0-t4.*(t6+t7-t18+t76+t77+t78+t79-t3.*t16.*theta_dot1)+t3.*(t9+t10+t11-t16.*(t27+t28)-t14.*t43+t13.*t54-t17.*t50+t15.*t58-t4.*t16.*theta_dot1)-(t3.*theta_dot1.*theta_dot2)./2.0;t65+t70+t75-t76-t77-t78-t79+t80-t3.*t5.*(9.81e2./2.0e2)-(t8.*t14)./2.0-(t8.*t17)./2.0-(t12.*t13)./2.0-(t12.*t15)./2.0-(t32.*theta_dot3)./2.0-(t36.*theta_dot4)./2.0-(t3.*theta_dot1.*theta_dot2)./2.0;t70-(t8.*t14)./2.0-(t12.*t13)./2.0-(t32.*theta_dot3)./2.0;t75-(t8.*t17)./2.0-(t12.*t15)./2.0-(t36.*theta_dot4)./2.0];
