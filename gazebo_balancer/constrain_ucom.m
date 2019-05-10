function u_com = constrain_ucom(u_com_des, params)
    % Use the CWC criterion to constrain the input to the feedback-linearized
    % model such that contact constraints are enforced.

    mu = 0.3;  % coefficient of friction

    % CasADi formulation
    import casadi.*
    
    % Optimization Variables
    f_c1 = MX.sym('f_c1',2,1);
    f_c2 = MX.sym('f_c2',2,1);
    u_com = MX.sym('u_com',2,1);

    % Equality (dynamic) constraints
    A_eq = params.Jcom*inv(params.H)*(params.Jcom'*[eye(2) eye(2)]);
    A_eq = [A_eq, -eye(2)];
    b_eq = -params.Jcom*inv(params.H)*(-params.C)-params.Jcom_dot*params.qdot;

    % Inequality (friction cone) constraints
    A_ineq = [ 1 -mu  0   0;
              -1 -mu  0   0;
               0   0  1 -mu;
               0   0 -1 -mu];
    A_ineq = [A_ineq, zeros(4,2)];

    A = [A_eq;A_ineq];               % Total constraints are a combination of these
    ubg = [b_eq; zeros(4,1)]; 
    lbg = [b_eq; -Inf(4,1)];

    qp = struct;
    qp.x = vertcat(f_c1,f_c2, u_com);            % decision variables
    qp.f = (u_com_des-u_com)'*(u_com_des-u_com); % objective function
    qp.g = A*[f_c1;f_c2;u_com];                        % constraints

    S = qpsol('S','qpoases',qp);
    res = S('lbg',lbg,'ubg',ubg);
    
    u_com = res.x(5:6);

end
