function [res,f_c1,f_c2] = check_force(f_com, com_Xf_c1, com_Xf_c2)
    % This function checks if the following force f_com in the x-y plane
    % is valid, returning 1 if yes and otherwise 0. Depends on the current
    % spatial force transforms to the center of mass 

    tic
    mu = 0.3;

    % CasADi formulation
    import casadi.*
    f_c1 = MX.sym('f_c1',2,1);
    f_c2 = MX.sym('f_c2',2,1);

    % Constraints
    A = [ eye(2) eye(2);    % Equality (sum of forces) constraint
          1 -mu  0   0;      % Inequality (cone) constraint
         -1 -mu  0   0;
          0   0  1 -mu;
          0   0 -1 -mu ];

    ubg = [f_com; zeros(4,1)];  % estabilish an equality and an inequality constraint
    lbg = [f_com; -Inf(4,1)];

    qp = struct;
    qp.x = vertcat(f_c1,f_c2);            % decision variables
    qp.f = 0;                                    % objective function
    qp.g = A*[f_c1;f_c2];                        % constraints


    try
        S = qpsol('S','qpoases',qp);
        res = S('lbg',lbg,'ubg',ubg);
        
        f_c1 = res.x(1:2);
        f_c2 = res.x(3:4);
        
        res = S.stats.success;
    catch e
        %disp(e)
        res = 0
    end

    % YALMIP formulation
    %f_c1 = sdpvar(2,1);
    %f_c2 = sdpvar(2,1);

    %F = [ [0;0;0;f_com;0] == com_Xf_c1*[0;0;0;f_c1;0] + com_Xf_c2*[0;0;0;f_c2;0] ];
    %F = F + [ f_c1(1) <= mu*f_c1(2) ];
    %F = F + [ f_c2(1) <= mu*f_c2(2) ];
    %F = F + [ f_c1(1) >= -mu*f_c1(2) ];
    %F = F + [ f_c2(1) >= -mu*f_c2(2) ];
    %
    %res = optimize(F,[]);
    %res = 1-res.problem;

    %f_c1 = value(f_c1)
    %f_c2 = value(f_c2)

    %% CVX formulation
    %cvx_begin quiet
    %    variable f_c1(2)  % x and y forces on contacts
    %    variable f_c2(2)

    %    minimize(0)
    %    subject to
    %        % forces must sum up to the net result
    %        %[0;0;0;f_com;0] == com_Xf_c1*[0;0;0;f_c1;0] + com_Xf_c2*[0;0;0;f_c2;0]
    %        f_com == f_c1 + f_c2;

    %        % forces must obey friction
    %        
    %        f_c1(1) <= mu*f_c1(2)
    %        f_c2(1) <= mu*f_c2(2)
    %        f_c1(1) >= -mu*f_c1(2)
    %        f_c2(1) >= -mu*f_c2(2)
    %        %norm(f_c1(1)) <= mu*f_c1(2)
    %        %norm(f_c2(1)) <= mu*f_c2(2)
    %cvx_end

    %f_c1
    %f_c2

    %if cvx_status=="Solved"
    %    res = 1;
    %else
    %    res = 0;
    %end

    toc
end
