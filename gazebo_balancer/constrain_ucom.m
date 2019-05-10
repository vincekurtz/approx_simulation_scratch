function u_com = constrain_ucom(u_com_des, params)
    % Use the CWC criterion to constrain the input to the feedback-linearized
    % model such that contact constraints are enforced.

    tic;
    mu = 0.3;  % coefficient of friction

    % CasADi formulation
    import casadi.*
    
    opti = casadi.Opti();
    u_com = opti.variable(2,1);
    f_c1 = opti.variable(2,1);
    f_c2 = opti.variable(2,1);

    % Equality (dynamic) constraints
    A_eq = params.Jcom*inv(params.H)*(params.Jcom'*[eye(2) eye(2)]);
    b_eq = params.Jcom*inv(params.H)*(-params.C)+params.Jcom_dot*params.qdot;

    % Inequality (friction cone) constraints
    A_ineq = [ 1 -mu  0   0;
              -1 -mu  0   0;
               0   0  1 -mu;
               0   0 -1 -mu ];

    opti.minimize((u_com-u_com_des)'*(u_com-u_com_des));

    opti.subject_to( A_eq*[f_c1;f_c2]+b_eq == u_com );
    opti.subject_to( A_ineq*[f_c1;f_c2] <= 0 );

    opti.solver('ipopt');
    sol = opti.solve();
    u_com = sol.value(u_com);
    
    %u_com = MX.sym('u_com',2,1);
    %f_c1 = MX.sym('f_c1',2,1);
    %f_c2 = MX.sym('f_c2',2,1);

    %% Equality (dynamic) constraints
    %A_eq = params.Jcom*inv(params.H)*(params.Jcom'*[eye(2) eye(2)]);
    %b_eq = params.Jcom*inv(params.H)*(-params.C)+params.Jcom_dot*params.qdot;

    %% Inequality (friction cone) constraints
    %A_ineq = [ 1 -mu  0   0;
    %          -1 -mu  0   0;
    %           0   0  1 -mu;
    %           0   0 -1 -mu ];

    %A = [A_eq;A_ineq];               % Total constraints are a combination of these
    %ubg = [u_com_des-b_eq; zeros(4,1)]; 
    %lbg = [u_com_des-b_eq; -Inf(4,1)];

    %qp = struct;
    %qp.x = vertcat(f_c1,f_c2);            % decision variables
    %qp.f = 0;%(u_com_des-u_com)'*(u_com_des-u_com); % objective function
    %qp.g = A*[f_c1;f_c2];                        % constraints

    %S = qpsol('S','qpoases',qp);
    %res = S('lbg',lbg,'ubg',ubg)
    %
    %u_com = res.x(1:2);

    toc


end
