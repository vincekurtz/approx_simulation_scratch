% Script to verify vector projection stuff

o_p_com = [0.6;0.5];  % (x,y) position of center of mass in O fame

o_p_c1 = [0.5;0];  % position of contacts in O frame
o_p_c2 = [-0.5;0];

% position of com in contact frames
c1_p_com = o_p_com-o_p_c1;
c2_p_com = o_p_com-o_p_c2;

f1 = [0;0.5];  % force 1
%f1_com = vector_proj(f1,c1_p_com);      % force 1 projected to the CoM

% Check if f1 is valid
[res,f_c1,f_c2] = check_force(f1-[0;2*-9.81], o_p_com)

figure;
hold on
scatter(o_p_com(1), o_p_com(2),'filled')
scatter(o_p_c1(1), o_p_c1(2),'filled')
scatter(o_p_c2(1), o_p_c2(2),'filled')

%plot_vector(o_p_c2, c2_p_com)
%plot_vector(o_p_c1, c1_p_com)

plot_vector(o_p_com,f1)
plot_vector(o_p_c1,f_c1)
plot_vector(o_p_c2,f_c2)

f_c1_proj = vector_proj(f_c1,c1_p_com);
f_c2_proj = vector_proj(f_c2,c2_p_com);


xlim([-1,1])
ylim([0,2])

function ahat = vector_proj(a,b)
    % project vector a onto vector b
    P = (b*b') / (b'*b);
    ahat = P*a;
    %ahat = dot(a,b)/norm(b)*b;
end

function plot_vector(frame, vector)
    % plot a vector from start_pos to end_pos
    x0 = frame(1);
    y0 = frame(2);

    x1 = x0 + vector(1);
    y1 = y0 + vector(2);

    plot([x0 x1], [y0 y1])
end
