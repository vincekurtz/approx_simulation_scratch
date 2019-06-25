function zmp = compute_ZMP(u_com, q, m, g)
    
    % Spatial momentum derivative in the CoM frame
    hd_G = [0;0;u_com;0];

    % Spatial force in the CoM frame
    f_G = hd_G - [0;0;0;0;-m*g;0];

    % Spatial force in the 0 frame
    f_0 = inv(Xf_0(q))*f_G;

    % ZMP position
    zmp = f_0(3)/f_0(5);
end
