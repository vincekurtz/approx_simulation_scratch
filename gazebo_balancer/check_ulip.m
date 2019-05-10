function res = check_ulip(u_lip, params)
    % Use the CWC criterion to check if the given input to the lip model
    % will preserve ground contacts. 

    mu = 0.3;

    % Compute the resulting force on the COM
    xdd_com = params.Rif*u_lip + params.Qif*params.x_lip + params.K*(params.x_com-params.x_lip);
    tau = InterfaceFcn(u_lip, params.x_lip, params.x);
    f_com = inv(params.Jcom')*tau; 

    % CasADi formulation
    import casadi.*
    f_c1 = MX.sym('f_c1',2,1);
    f_c2 = MX.sym('f_c2',2,1);

    % Constraints
    Rbar = 1/params.Rif;
    u_lip
    Rbar*(xdd_com - params.Qif*params.x_lip - params.K*(params.x_com-params.x_lip))

    A_eq = params.Jcom*inv(params.H)*(params.Jcom'*[eye(2) eye(2)]);
    b_eq = params.Jcom*inv(params.H)*(-params.C)+params.Jcom_dot*params.qdot;

    A_ineq = [ 1 -mu  0   0;      % Inequality (cone) constraint
              -1 -mu  0   0;
               0   0  1 -mu;
               0   0 -1 -mu ];

    A = [A_eq;A_ineq];               % Total constraints are a combination of these
    ubg = [xdd_com-b_eq; zeros(4,1)]; 
    lbg = [xdd_com-b_eq; -Inf(4,1)];

    qp = struct;
    qp.x = vertcat(f_c1,f_c2);      % decision variables
    qp.f = 0;                       % objective function
    qp.g = A*[f_c1;f_c2];           % constraints


    try
        S = qpsol('S','qpoases',qp);
        res = S('lbg',lbg,'ubg',ubg)
        
        f_c1 = res.x(1:2);
        f_c2 = res.x(3:4);

        res = S.stats.success;
    catch e
        %disp(e)
        res = 0;
    end

end
