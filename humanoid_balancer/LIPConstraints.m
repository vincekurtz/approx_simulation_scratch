function cons = LIPConstraints(params)
    % Returns an object which respresents linear constraints and a quadratic 
    % objective for MPC with the LIP model. These constrains enforce both contacts
    % and the simulation relation with the full-order model.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cost Function                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Running Cost
    Q = diag([10;0;0;10;0]);
    R = 5.00;

    % Terminal cost
    Qf = 100*Q;

    % Running and terminal costs together
    QQ = kron(eye(params.N-1),Q);  % matrix with Q on diagonal

    % Cost function is given by x_lip(:)'QQ*x_lip(:) + u_lip(:)'*RR*u_lip(:)
    cons.QQ = blkdiag(QQ,Qf);           % the last element is the terminal cost
    cons.RR = kron(eye(params.N-1),R);  % matrix with R on diagonal

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial Conditions Constraint    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A_init_com = [eye(5) zeros(5, 5*(params.N-1) + 3*(params.N-1) + 5*params.N + 1*(params.N-1))];
    A_init_lip = [zeros(5, 5*params.N + 3*(params.N-1)) eye(5) zeros(5, 5*(params.N-1) + 1*(params.N-1))];

    % Initial conditions constraint is A_init*[x_com(:);u_com(:);x_lip(:);u_lip(:)] == [x_com_init;u_com_init]
    cons.A_init = [A_init_com;
                   A_init_lip];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dynamics Constraint              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % Dynamic constraint is A_dynamics*[x_com(:);u_com(:);x_lip(:);u_lip(:)] == b_dynamics
    cons.A_dynamics = [A_com_dynamics                  zeros(5*(params.N-1), 5*params.N+1*(params.N-1));
                       zeros(5*(params.N-1), 5*params.N+3*(params.N-1))                  A_lip_dynamics];
    cons.b_dynamics = zeros(5*(params.N-1)+5*(params.N-1),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Constraint             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Interface constraint is A_interface*[x_com(:);u_com(:);x_lip(:);u_lip(:)] == b_interface
    cons.A_interface = [kron(eye(params.N-1), params.K)           ...
                        zeros(3*(params.N-1),5)                   ...
                        kron(eye(params.N-1), -eye(3))            ...
                        kron(eye(params.N-1), params.Q-params.K)  ...
                        zeros(3*(params.N-1),5)                   ...
                        kron(eye(params.N-1), params.R)];
    cons.b_interface = zeros(3*(params.N-1),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Accelration Bound Constraint     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bnd = 5.0;
    A_ucom_bound = [ [eye(2);-eye(2)] zeros(4,1) ];
    A_ucom_bound = [kron(eye(params.N-1), A_ucom_bound) ];
    
    % Accelration bound constraint is A_bnd*[x_com(:);u_com(:);x_lip(:);u_lip(:)] <= b_bnd
    cons.A_bnd = [ zeros(4*(params.N-1), 5*params.N), A_ucom_bound, zeros(4*(params.N-1), 5*(params.N) + 1*(params.N-1))];
    cons.b_bnd = bnd*ones(4*(params.N-1),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contact Constraints              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    % Contact constraint is A_contact*[x_com(:);u_com(:);x_lip(:);u_lip(:)] <= b_contact
    cons.A_contact = [A_contact_1;
                      A_contact_2;
                      A_contact_3;
                      A_contact_4];
    cons.b_contact = repmat(b,4*(params.N-1),1);

    % Include the number of timesteps
    cons.N = params.N;

end

function y = S(x)
    % Skew symmetric cross product matrix
    y = [0     -x(3)   x(2);
         x(3)   0     -x(1);
         -x(2)  x(1)    0  ];
end
