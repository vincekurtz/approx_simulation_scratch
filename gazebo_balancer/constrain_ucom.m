function u_com = constrain_ucom(u_com_des, params)
    % Use the CWC criterion to constrain the input to the feedback-linearized
    % model such that contact constraints are enforced.

    import casadi.*
    
    % Optimization Variables
    f_c1 = MX.sym('f_c1',2,1);
    f_c2 = MX.sym('f_c2',2,1);
    u_com = MX.sym('u_com',2,1);

    % Equality (dynamic) constraints
    A_eq = params.Jcom*inv(params.H)*(params.Jcom'*[eye(2) eye(2)]);  % sum of forces constraint
    A_eq = [A_eq, -eye(2)]; 
    b_eq = params.Jcom*inv(params.H)*(params.C)-params.Jcom_dot*params.qdot;

    S = [0 0 1 0 0 0];      % "cross product" constraint (resulting torque on COM must be zero)
    C = [zeros(3,2);eye(2);zeros(1,2)];
    A_cross = [S*params.X1*C S*params.X2*C zeros(1,2)];
    b_cross = 0;
    A_eq = [A_eq; A_cross];
    b_eq = [b_eq; 0];

    % Inequality (friction cone) constraints
    mu = params.mu;  % coefficient of friction
    A_ineq = [ 1 -mu  0   0;
              -1 -mu  0   0;
               0   0  1 -mu;
               0   0 -1 -mu];
    A_ineq = [A_ineq, zeros(4,2)];

    A = [A_eq;A_ineq];               % Total constraints are a combination of these
    ubg = [b_eq; zeros(4,1)]; 
    lbg = [b_eq; -Inf(4,1)];

    % QP setup
    qp = struct;
    qp.x = vertcat(f_c1,f_c2, u_com);            % decision variables

    % Objective function minimizes simulation function V(x)
    x_lip_plus = params.x_lip + (params.A_lip*params.x_lip + params.B_lip*params.u_lip)*params.dt;
    x_com_plus = params.x_com + (params.A_com*params.x_com + params.B_com*u_com)*params.dt;
    qp.f = (x_com_plus-x_lip_plus)'*params.M*(x_com_plus-x_lip_plus);

    % Alternative objective function is to minimize the distance to the nominal control.
    qp.f = (u_com_des - u_com)'*(u_com_des - u_com);

    % Constraints to ensure that we remain in contact
    qp.g = A*[f_c1;f_c2;u_com];

    % Solve the QP
    S = qpsol('S','qpoases',qp);
    res = S('lbg',lbg,'ubg',ubg);
    
    u_com = full(res.x(5:6));  % convert to double


end
