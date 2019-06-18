function u_lip = GenerateLIPTrajectory(x_lip_init, x_com_init, params)
    % Compute a trajectory for the LIP model that respects contact constraints
    % of the full model.

    import casadi.*;
   
    % Optimization variables
    x_com = SX.sym('x_com',5*params.N);
    u_com = SX.sym('u_com',3*(params.N-1));
    x_lip = SX.sym('x_lip',5*params.N);
    u_lip = SX.sym('u_lip',1*(params.N-1));
    qp.x = vertcat( x_com, u_com, x_lip, u_lip );

    % Cost Function
    qp.f = x_lip'*params.QQ*x_lip + u_lip'*params.RR*u_lip;

    % Initial Condition Constraints depend on supplied initial conditions
    b_init = [x_com_init;
              x_lip_init];

    % Equality constraints from dynamics, initial conditions, and the interface
    A_eq = [params.A_init;
            params.A_dynamics;
            params.A_interface];
    b_eq = [b_init;
            params.b_dynamics;
            params.b_interface];
    
    % Inequality constraints from acceleration bound and contact constraints
    A_ineq = [params.A_bnd;
              params.A_contact];  
    b_ineq = [params.b_bnd;
              params.b_contact];

    % Putting constraints in standard form
    qp.g = [A_eq;A_ineq]*[x_com;u_com;x_lip;u_lip];
    lbg = [b_eq;-Inf(size(b_ineq))];
    ubg = [b_eq;b_ineq];

    % Solve the QP
    options.ipopt.print_level = 0;
    options.print_time = false;
    Solver = nlpsol('Solver','ipopt',qp, options);
    %Solver = qpsol('S','qpoases',qp);
    r = Solver('lbg',lbg,'ubg',ubg);
    x_opt = full(r.x);
    u_lip = x_opt(end-params.N+2:end)';

end

function y = S(x)
    % Skew symmetric cross product matrix
    y = [0     -x(3)   x(2);
         x(3)   0     -x(1);
         -x(2)  x(1)    0  ];
end
