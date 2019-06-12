function u_lip = LIPController(x_lip, x_com, params)
    % Compute a control for the (simplified) LIP model
   
    import casadi.*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameterizing the Cost Function %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % LIP dynamics
    A_lip = params.A_lip;
    B_lip = params.B_lip;

    % Running Cost
    Q = diag([100;10;0;0;0]);
    R = 0.01;

    % Terminal cost based on optimal cost-to-go of LQR
    %[K,Qf] = lqr(A_lip,B_lip,Q,R)
    Qf = [31.8676 0.0772 0 0 0;
           0.0772 0.0244 0 0 0;
           0      0      0 0 0;
           0      0      0 0 0;
           0      0      0 0 0];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameterizing Contact Constraints  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Force transform from the CoM frame to the 0 frame
    o_Xf_com = [eye(3)  [0         0        x_com(2) ;
                         0         0        -x_com(1);
                         -x_com(2) x_com(1) 0       ];
                zeros(3)   eye(3)                     ];

    % Constraints on the contact wrench (note that we only consider motion in the x-y plane)
    mu = params.friction_coef;
    l = params.foot_width;
    A = [0 0  0  0  -1 0;   % positive normal force
         0 0  0  1 -mu 0;   % Coulomb friction
         0 0  0 -1 -mu 0;
         0 0  1  0  -l 0;   % Center of pressure constraint
         0 0 -1  0  -l 0];

    % Constraints on u_com
    A_cwc = A*o_Xf_com;
    A_cwc = A_cwc*[zeros(2,3);eye(3);zeros(1,3)]; % map u_com to [0;0;u_com;0] implicitly
    b_cwc = A*o_Xf_com*[0;0;0;0;params.mg;0];

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization via Direct Collocation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K_lip = [ -11 -2.1 ];
    u_lip_des = -K_lip*[x_lip(1);1/params.m*x_lip(4)];

    opti = casadi.Opti();
    u_lip = opti.variable(1);
    u_com = opti.variable(3,1);

    opti.minimize((u_lip - u_lip_des)'*(u_lip - u_lip_des));
    opti.subject_to(u_com == params.R*u_lip + params.Q*x_lip + params.K*(x_com- x_lip));
    opti.subject_to( A_cwc*u_com <= b_cwc );

    options.ipopt.print_level = 1;
    options.print_time = true;
    opti.solver('ipopt', options);
    sol = opti.solve();
    u_lip = sol.value(u_lip(1));


end

