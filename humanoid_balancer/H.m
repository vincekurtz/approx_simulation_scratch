function H_sym = H(in1)
%H
%    H_SYM = H(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    18-Jun-2019 14:53:04

theta2 = in1(2,:);
theta3 = in1(3,:);
theta4 = in1(4,:);
t2 = cos(theta4);
t3 = t2./2.0;
t4 = sin(theta4);
t5 = cos(theta3);
t6 = t2.^2;
t7 = t4.^2;
t8 = sin(theta3);
t9 = t6+t7+2.0;
t10 = t2+1.0./2.0;
t11 = sin(theta2);
t12 = t5.^2;
t13 = t9.*t12;
t14 = t8.^2;
t15 = t9.*t14;
t16 = t3+t6+t7+2.0;
t17 = t5.*t16;
t18 = t5.*t9;
t19 = t2.*t10;
t20 = t7+t18+t19+2.0;
t21 = t2.*t4;
t22 = t8.*t9;
t26 = t4.*t10;
t23 = t21+t22-t26;
t24 = cos(theta2);
t25 = t5.*t20;
t27 = t8.*t23;
t28 = t13+t15+1.0;
t34 = (t4.*t8)./2.0;
t29 = t13+t15+t17-t34+1.0./2.0;
t30 = t24.*t29;
t31 = (t4.*t5)./2.0;
t32 = t8.*t16;
t33 = t31+t32;
t35 = (t2.*t5)./2.0;
t38 = t11.*t33;
t36 = t3+t7+t17+t19+t25+t27+t30-t34-t38+2.753766666666667;
t37 = t24.*(t17-t34);
t39 = t3+t7+t17+t19-t34+2.419183333333333;
t40 = t34-t35;
t41 = (t2.*t8)./2.0;
t42 = t31+t41;
t43 = t3-t34+t35-t11.*t42-t24.*t40+3.345833333333333e-1;
t44 = t3-t34+t35+3.345833333333333e-1;
t45 = t3+3.345833333333333e-1;
H_sym = reshape([t3+t7+t17+t19+t25+t27+t30+t11.*(t5.*t23-t8.*t20+t11.*t28)-(t4.*t8)./2.0-t11.*t33+t24.*(t25+t27+t24.*t28+1.0./2.0)+3.08835,t36,t3+t7+t17+t19-t34+t37-t38+2.419183333333333,t43,t36,t3+t7+t17+t19+t25+t27-t34+2.753766666666667,t39,t44,t3+t7+t17+t19-t34+t37-t11.*t33+2.419183333333333,t39,t3+t7+t19+2.419183333333333,t45,t43,t44,t45,3.345833333333333e-1],[4,4]);
