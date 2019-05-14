function [res] = check_force2(f_com, X1, X2)
    % This function checks if the following force f_com in the x-y plane
    % is valid, returning 1 if yes and otherwise 0

    mu = 0.3;

    % CasADi formulation
    import casadi.*

    opti = casadi.Opti();
    f_c1 = opti.variable(2,1);
    f_c2 = opti.variable(2,1);

    % Force constraint
    S = [0 0 1 0 0 0;   % Consider only the n_z, f_x, and f_y constraints since
         0 0 0 1 0 0;   % we're restricted to the x-y plane
         0 0 0 0 1 0];
    opti.subject_to( S*X1*[0;0;0;f_c1;0] + S*X2*[0;0;0;f_c2;0] == [0;f_com] );
    %opti.subject_to( f_c1 + f_c2 == [f_com] );

    % Friction constraints
    A = [1 -mu ; -1 -mu];
    opti.subject_to( A*f_c1 <= 0 );
    opti.subject_to( A*f_c2 <= 0);

    try
        opti.solver('ipopt');
        sol = opti.solve();

        f_c1 = sol.value(f_c1);
        f_c2 = sol.value(f_c2);
        f_com
        X1*[0;0;0;f_c1;0] + X2*[0;0;0;f_c2;0]

        res = sol.isvalid;
    catch e
        disp(e)
        res = 0;
    end


end
