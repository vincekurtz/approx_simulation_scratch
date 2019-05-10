function [res,f_c1,f_c2] = check_force(f_com)
    % This function checks if the following force f_com in the x-y plane
    % is valid, returning 1 if yes and otherwise 0

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

end
