function Ad_com_qd_sym = Ad_com_qd(in1,in2)
%AD_COM_QD
%    AD_COM_QD_SYM = AD_COM_QD(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Jun-2019 10:22:15

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
t51 = t22.*theta_dot3;
t52 = t12.*t14;
t43 = t37+t39+t42-t51-t52;
t44 = t8.*t15;
t46 = t23.*t26;
t47 = theta_dot4./2.0;
t48 = t27+t28+t47;
t49 = t23.*t48;
t53 = t26.*theta_dot4;
t54 = t12.*t17;
t50 = t44+t46+t49-t53-t54;
t55 = t27+t28;
t56 = t8.*t14;
t57 = t12.*t13;
t58 = t32.*theta_dot3;
t67 = t19.*t32;
t59 = t56+t57+t58-t67;
t60 = t13.*t59;
t61 = t8.*t17;
t62 = t12.*t15;
t63 = t36.*theta_dot4;
t68 = t23.*t36;
t64 = t61+t62+t63-t68;
t65 = t15.*t64;
t69 = t14.*t43;
t70 = t17.*t50;
t71 = t16.*t55;
t72 = t4.*t16.*theta_dot1;
t66 = t9+t10+t11+t60+t65-t69-t70-t71-t72;
t73 = t3.*t66;
t74 = t13.*t43;
t75 = t15.*t50;
t76 = t14.*t59;
t77 = t17.*t64;
t79 = t3.*t16.*theta_dot1;
t78 = t6+t7-t18+t74+t75+t76+t77-t79;
t80 = t2.*(9.81e2./1.0e2);
t90 = t4.*t78;
t81 = t73+t80-t90;
t82 = t5.*(9.81e2./1.0e2);
t83 = t4.*t66;
t84 = theta_dot1.^2;
t85 = t3.*t78;
t92 = t84./2.0;
t86 = t82+t83+t85-t92;
t87 = theta1+theta2+theta3;
t88 = theta1+theta2+theta4;
t89 = theta1+theta2;
t91 = t5.*t81;
t93 = t5.*t86;
t94 = t2.*t81;
Ad_com_qd_sym = [t2.*(9.81e2./2.0e2)+t73-t76-t77+t2.*t4.*(9.81e2./2.0e2)-t3.*t5.*(9.81e2./2.0e2)-(t8.*t14)./2.0-(t8.*t17)./2.0-(t12.*t13)./2.0-(t12.*t15)./2.0-t13.*t43-t15.*t50-t4.*t78-(t32.*theta_dot3)./2.0-(t36.*theta_dot4)./2.0-(t93+t94).*(t2.*(7.0./8.0)+cos(t87)./8.0+cos(t88)./8.0+cos(t89).*(5.0./8.0))-(t91-t2.*t86).*(t5.*(7.0./8.0)+sin(t87)./8.0+sin(t88)./8.0+sin(t89).*(5.0./8.0))+t19.*((t13.*t16)./2.0+(t3.*t14.*theta_dot1)./2.0+(t4.*t13.*theta_dot1)./2.0)+t23.*((t15.*t16)./2.0+(t4.*t15.*theta_dot1)./2.0+(t3.*t17.*theta_dot1)./2.0)+(t3.*t16.*theta_dot1)./2.0-(t3.*theta_dot1.*theta_dot2)./2.0;-t91+t2.*t86;t93+t94-9.81e2./2.5e1];
