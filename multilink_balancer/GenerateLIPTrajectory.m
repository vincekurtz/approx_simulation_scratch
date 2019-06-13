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

    opti = casadi.Opti();
    
    u_lip = opti.variable(1,params.N-1);
    x_lip = opti.variable(5,params.N);
    u_com = opti.variable(3,params.N-1);
    x_com = opti.variable(5,params.N);

    % Cost Function
    QQ = kron(eye(params.N-1),Q);  % matrix with Q on diagonal
    QQ = blkdiag(QQ,Qf);           % the last element is the terminal cost
    RR = kron(eye(params.N-1),R);  % matrix with R on diagonal
    J = x_lip(:)'*QQ*x_lip(:) + u_lip(:)'*RR*u_lip(:);

    % Constraints
    opti.subject_to(x_lip(:,1) == x_lip_init);
    opti.subject_to(x_com(:,1) == x_com_init);

    % Constraints on the contact wrench (note that we only consider motion in the x-y plane)
    mu = 0.2;
    l = 0.5;
    mg = -39.24;
    A = [0 0  0  0  -1 0;   % positive normal force
         0 0  0  1 -mu 0;   % Coulomb friction
         0 0  0 -1 -mu 0;
         0 0  1  0  -l 0;   % Center of pressure constraint
         0 0 -1  0  -l 0];
   
    % Interface Constraint
    A_interface_eq = [kron(eye(params.N-1), params.R)           ...
                      kron(eye(params.N-1), params.Q-params.K)  ...
                      kron(eye(params.N-1), -eye(3))            ...
                      kron(eye(params.N-1), params.K)];
    opti.subject_to( A_interface_eq*[u_lip(:); x_lip(1:5*(params.N-1))'; u_com(:); x_com(1:5*(params.N-1))'] == 0 );

    % Acceleration bound (for linearizing contact constraints)
    bnd = 3.0;
    opti.subject_to( -bnd <= u_com(2,:) <= bnd );  % bound linear acceleration
    opti.subject_to( -bnd <= u_com(3,:) <= bnd );

    % Contact constraints
    Au = A*[zeros(2,3);eye(3);zeros(1,3)];
    Ax_1 = A*[S([0;mg;0])-S([bnd;bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_2 = A*[S([0;mg;0])-S([-bnd;bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_3 = A*[S([0;mg;0])-S([bnd;-bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    Ax_4 = A*[S([0;mg;0])-S([-bnd;-bnd;0]);zeros(3)]*[eye(2) zeros(2,3); zeros(1,5)];
    b = A*[0;0;0;0;mg;0];

    A_contact_1 = [kron(eye(params.N-1),Au) kron(eye(params.N-1),Ax_1)];
    A_contact_2 = [kron(eye(params.N-1),Au) kron(eye(params.N-1),Ax_2)];
    A_contact_3 = [kron(eye(params.N-1),Au) kron(eye(params.N-1),Ax_3)];
    A_contact_4 = [kron(eye(params.N-1),Au) kron(eye(params.N-1),Ax_4)];
    b_contact = repmat(b,params.N-1,1);

    opti.subject_to( A_contact_1*[u_com(:);x_com(1:5*(params.N-1))'] <= b_contact )
    opti.subject_to( A_contact_2*[u_com(:);x_com(1:5*(params.N-1))'] <= b_contact )
    opti.subject_to( A_contact_3*[u_com(:);x_com(1:5*(params.N-1))'] <= b_contact )
    opti.subject_to( A_contact_4*[u_com(:);x_com(1:5*(params.N-1))'] <= b_contact )

    % Dynamic (forward Euler) constraints
    A_lip_eq = [eye(5)-params.A_lip*params.dt -eye(5) zeros(5,5*(params.N-2))];
    B_lip_eq = [params.B_lip*params.dt zeros(5,params.N-2)];

    A_com_eq = [eye(5)-params.A_com*params.dt -eye(5) zeros(5,5*(params.N-2))];
    B_com_eq = [params.B_com*params.dt zeros(5,3*(params.N-2))];

    for t=2:params.N-1
        A_lip_last_row = A_lip_eq((end-4):end,:);
        A_lip_eq = [A_lip_last_row;
                    circshift(A_lip_last_row,[0,5])];
        B_lip_last_row = B_lip_eq((end-4):end,:);
        B_lip_eq = [B_lip_last_row;
                    circshift(B_lip_last_row,[0,1])];

        A_com_last_row = A_com_eq((end-4):end,:);
        A_com_eq = [A_com_last_row;
                    circshift(A_com_last_row,[0,5])];
        B_com_last_row = B_com_eq((end-4):end,:);
        B_com_eq = [B_com_last_row;
                    circshift(B_com_last_row,[0,3])];
    end

    opti.subject_to( [A_lip_eq B_lip_eq]*[x_lip(:);u_lip(:)] == 0 )
    opti.subject_to( [A_com_eq B_com_eq]*[x_com(:);u_com(:)] == 0 )

    % Run the optimization
    opti.minimize(J);

    options.ipopt.print_level = 0;
    options.print_time = false;
    opti.solver('ipopt', options);
    %opti.solver('sqpmethod');
    sol = opti.solve();
    u_lip = sol.value(u_lip);

end

function y = S(x)
    % Skew symmetric cross product matrix
    y = [0     -x(3)   x(2);
         x(3)   0     -x(1);
         -x(2)  x(1)    0  ];
end
