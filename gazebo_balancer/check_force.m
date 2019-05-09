function [res,f_c1,f_c2] = check_force(f_com, p_com)
    % This function checks if the following force f_com in the x-y plane
    % is valid, returning true if yes and otherwise no. Depends on the current
    % position of the center of mass, p_com

    p_c1 = [0.5;0];  % positions of foot contacts
    p_c2 = [-0.5;0]; 

    % position of com in the contact frames
    c1_p_com = p_com - p_c1;
    c2_p_com = p_com - p_c2;

    % Projection matrix onto the vector from the contact to the center of mass
    P1 = (c1_p_com*c1_p_com')/(c1_p_com'*c1_p_com);
    P2 = (c2_p_com*c2_p_com')/(c2_p_com'*c2_p_com);

    mu = 0.3;
    cvx_begin quiet
        variable f_c1(2)  % x and y forces on contacts
        variable f_c2(2)

        minimize(0)
        subject to
            % forces must sum up to the net result
            f_com == P1*f_c1 + P2*f_c2;

            % forces must obey friction
            norm(f_c1(1)) <= mu*f_c1(2)
            norm(f_c2(1)) <= mu*f_c2(2)
    cvx_end

    if cvx_status=="Solved"
        f_c1;
        f_c2;
        f_com = P1*f_c1 + P2*f_c2;
        res = 1;
    else
        res = 0;
    end
end
