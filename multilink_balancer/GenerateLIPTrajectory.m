function u_lip = GenerateLIPTrajectory(x_lip_init, x_com_init, params)
    % Compute a trajectory for the LIP model that respects contact constraints
    % of the full model.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameterizing the Cost Function %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Running Cost
    Q = diag([100;0;0;10;0]);
    R = 0.01;

    % Terminal cost based on optimal cost-to-go of LQR
    Qf = [32.1091 0 0  0.1544 0;
           0      0 0  0      0;
           0      0 0  0      0;
           0.1544 0 0  0.0491 0;
           0      0 0  0      0];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization via Direct Collocation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import casadi.*;
   
    % Optimization variables
    x_com = SX.sym('x_com',5,params.N);
    u_com = SX.sym('u_com',3,params.N-1);
    x_lip = SX.sym('x_lip',5,params.N);
    u_lip = SX.sym('u_lip',1,params.N-1);
    qp.x = vertcat( x_com(:), u_com(:), x_lip(:), u_lip(:) );

    % Cost Function
    QQ = kron(eye(params.N-1),Q);  % matrix with Q on diagonal
    QQ = blkdiag(QQ,Qf);           % the last element is the terminal cost
    RR = kron(eye(params.N-1),R);  % matrix with R on diagonal
    qp.f = x_lip(:)'*QQ*x_lip(:) + u_lip(:)'*RR*u_lip(:);

    % Initial Condition Constraints
    A_init_com = [eye(5) zeros(5, 5*(params.N-1) + 3*(params.N-1) + 5*params.N + 1*(params.N-1))];
    A_init_lip = [zeros(5, 5*params.N + 3*(params.N-1)) eye(5) zeros(5, 5*(params.N-1) + 1*(params.N-1))];
    A_init = [A_init_com;
              A_init_lip];
    b_init = [x_com_init;
              x_lip_init];

    % Dynamic (forward Euler) Constraints
    A_lip_eq = [eye(5)+params.A_lip*params.dt -eye(5) zeros(5,5*(params.N-2))];
    B_lip_eq = [params.B_lip*params.dt zeros(5,params.N-2)];

    A_com_eq = [eye(5)+params.A_com*params.dt -eye(5) zeros(5,5*(params.N-2))];
    B_com_eq = [params.B_com*params.dt zeros(5,3*(params.N-2))];

    for t=2:params.N-1
        A_lip_last_row = A_lip_eq((end-4):end,:);
        A_lip_eq = [A_lip_eq;
                    circshift(A_lip_last_row,[0,5])];

        B_lip_last_row = B_lip_eq((end-4):end,:);
        B_lip_eq = [B_lip_eq;
                    circshift(B_lip_last_row,[0,1])];

        A_com_last_row = A_com_eq((end-4):end,:);
        A_com_eq = [A_com_eq;
                    circshift(A_com_last_row,[0,5])];
        B_com_last_row = B_com_eq((end-4):end,:);
        B_com_eq = [B_com_eq;
                    circshift(B_com_last_row,[0,3])];
    end

    A_lip_dynamics = [A_lip_eq B_lip_eq];
    A_com_dynamics = [A_com_eq B_com_eq];
    A_dynamics = [A_com_dynamics                  zeros(5*(params.N-1), 5*params.N+1*(params.N-1));
                  zeros(5*(params.N-1), 5*params.N+3*(params.N-1))                  A_lip_dynamics];
    b_dynamics = zeros(5*(params.N-1)+5*(params.N-1),1);

    % Interface Constraints
    A_interface = [kron(eye(params.N-1), params.K)           ...
                   zeros(3*(params.N-1),5)                   ...
                   kron(eye(params.N-1), -eye(3))            ...
                   kron(eye(params.N-1), params.Q-params.K)  ...
                   zeros(3*(params.N-1),5)                   ...
                   kron(eye(params.N-1), params.R)];
    b_interface = zeros(3*(params.N-1),1);

    % Acceleration Bound Constraints
    bnd = 3.0;
    A_ucom_bound = [ [eye(2);-eye(2)] zeros(4,1) ];
    A_ucom_bound = [kron(eye(params.N-1), A_ucom_bound) ];
    A_bnd = [ zeros(4*(params.N-1), 5*params.N), A_ucom_bound, zeros(4*(params.N-1), 5*(params.N) + 1*(params.N-1))];
    b_bnd = bnd*ones(4*(params.N-1),1);

    % Contact Constraints
    mu = 0.2;
    l = 0.5;
    mg = -39.24;
    A = [0 0  0  0  -1 0;   % positive normal force
         0 0  0  1 -mu 0;   % Coulomb friction
         0 0  0 -1 -mu 0;
         0 0  1  0  -l 0;   % Center of pressure constraint
         0 0 -1  0  -l 0];
    
    Au = A*[zeros(2,3);eye(3);zeros(1,3)];
    Ax_1 = A*[S([0;mg;0])-S([bnd;bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_2 = A*[S([0;mg;0])-S([-bnd;bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_3 = A*[S([0;mg;0])-S([bnd;-bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_4 = A*[S([0;mg;0])-S([-bnd;-bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    b = A*[0;0;0;0;mg;0];

    A_contact_1 = [kron(eye(params.N-1),Ax_1) ...
                   zeros(5*(params.N-1),5)    ...
                   kron(eye(params.N-1),Au)   ...
                   zeros(5*(params.N-1),5*params.N+1*(params.N-1))];
    A_contact_2 = [kron(eye(params.N-1),Ax_2) ...
                   zeros(5*(params.N-1),5)    ...
                   kron(eye(params.N-1),Au)   ...
                   zeros(5*(params.N-1),5*params.N+1*(params.N-1))];
    A_contact_3 = [kron(eye(params.N-1),Ax_3) ...
                   zeros(5*(params.N-1),5)    ...
                   kron(eye(params.N-1),Au)   ...
                   zeros(5*(params.N-1),5*params.N+1*(params.N-1))];
    A_contact_4 = [kron(eye(params.N-1),Ax_4) ...
                   zeros(5*(params.N-1),5)    ...
                   kron(eye(params.N-1),Au)   ...
                   zeros(5*(params.N-1),5*params.N+1*(params.N-1))];

    A_contact = [A_contact_1;
                 A_contact_2;
                 A_contact_3;
                 A_contact_4];
    b_contact = repmat(b,4*(params.N-1),1);

    % Putting constraints together
    A_eq = [A_init;     % Equality constraints from dynamics, initial conditions, and the interface
            A_dynamics;
            A_interface];
    b_eq = [b_init;
            b_dynamics;
            b_interface];

    A_ineq = [A_bnd;        % Inequality constraints from acceleration bound and contact constraints
              A_contact];  
    b_ineq = [b_bnd;
              b_contact];

    qp.g = [A_eq;A_ineq]*[x_com(:);u_com(:);x_lip(:);u_lip(:)];


    lbg = [b_eq;-Inf(size(b_ineq))];
    ubg = [b_eq;b_ineq];

    % Solve the QP
    S = qpsol('S','qpoases',qp);
    %S = nlpsol('S','ipopt',qp);
    r = S('lbg',lbg,'ubg',ubg);
    x_opt = full(r.x);
    u_lip = x_opt(end-params.N+2:end)';

end

function y = S(x)
    % Skew symmetric cross product matrix
    y = [0     -x(3)   x(2);
         x(3)   0     -x(1);
         -x(2)  x(1)    0  ];
end
