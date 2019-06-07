function u_com = constrain_ucom(u_com_des, q, qd, params)
    % Use the CWC criterion to constrain the input to the feedback-linearized
    % model such that contact constraints are enforced.
    
    import casadi.*

    % Force transform from the CoM frame to the 0 frame
    o_Xf_com = inv(Xf_0(q));

    % Constraints on the contact wrench (note that we only consider motion in the x-y plane)
    mu = params.mu;
    l = 0.5;  % foot width
    A = [0 0  0  0  -1 0;   % positive normal force
         0 0  0  1 -mu 0;   % Coulomb friction
         0 0  0 -1 -mu 0;
         0 0  1  0  -l 0;   % Center of pressure constraint
         0 0 -1  0  -l 0];

    % Constraints on u_com
    A_cwc = A*o_Xf_com;
    A_cwc = A_cwc*[zeros(2,3);eye(3);zeros(1,3)]; % map u_com to [0;0;u_com;0] implicitly
    b_cwc = A*o_Xf_com*[0;0;0;0;params.mg;0];

    % Set up an optimization problem (QP) to constrain u_com
    u_com = MX.sym('u_com',3,1);

    qp = struct;
    qp.x = u_com;                                  % decision variables
    qp.f = (u_com-u_com_des)'*(u_com-u_com_des);  % objective

    qp.g = A_cwc*u_com;                            % constraints
    lbg = -Inf(5,1);
    ubg = b_cwc;

    % Solve the QP
    S = qpsol('S','qpoases',qp);
    res = S('lbg',lbg,'ubg',ubg);
    u_com = full(res.x);

end
