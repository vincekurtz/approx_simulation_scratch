function u_com = constrain_ucom(u_com_des, x_com)
    % Use the CWC criterion to constrain the input to the feedback-linearized
    % model such that contact constraints are enforced.
    
    import casadi.*

    % Set up an optimization problem (QP) to constrain u_com
    u_com = MX.sym('u_com',3,1);

    qp = struct;
    qp.x = u_com;                                  % decision variables
    qp.f = (u_com-u_com_des)'*(u_com-u_com_des);  % objective

    % Constraints
    qp.g = A_cwc(x_com)*u_com;  % note that A_cwc(x_com) = A_cwc(p_com) because x_com(1:2) = p_com
    lbg = -Inf(5,1);
    ubg = b_cwc(x_com);

    % Solve the QP
    S = qpsol('S','qpoases',qp);
    res = S('lbg',lbg,'ubg',ubg);
    u_com = full(res.x);

end
