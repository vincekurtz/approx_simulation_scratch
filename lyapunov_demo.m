%% Finding Local Lyapunov Certificats with SOS and YALMIP
% Adopted from Lemmon's Advanced Control Lecture Notes

% System: Van der Pool Oscillator
% x' = f(x)
% x = [x1 ; x2]
% x1' = -x2'
% x2' = x1-x2+(x1^2)x2

sdpvar x1 x2
x = [x1; x2];
f = [-x2 ;
     x1-x2+x1^2*x2];
 
% Goal: find a Lyapunov function V(x) 
% that locally certifies stability in the region X = {x | x1^2+x2^2 <= r}

% First we define the constraint g(x) such that X = {x | g(x) >= 0 }
r = 2.81;
g = r - (x1^2+x2^2);

% Now our goal is to find V(x) and s(x) such that
%   s(x) is SOS
%   V - eps(x'x) is SOS
%   -nablaV*f - s(x)g(x)-epsilon(x'x) is SOS
%
% If we can find such a V(x), s(x) for small epsilon, then V(x) is a 
% Lyapunov certificat for the region X, based on the positivestellensatz.

% First we define V(x) symbolically. We have to choose the structure 
% of it, which we do with monolist
v = monolist(x,4);  % this gives a vector of monomials up to degree d
n = length(v);
c = sdpvar(n,1);    % this is a vector of constants
V = c'*v;           % this is our Lyapunov function, with undefined constants

% Then we can similarly define nablaV, the gradient of V w.r.t x
nablaV = jacobian(V,x);
cs = sdpvar(n,1);
s = cs'*v;

% Here's our choice of epsilon
eps = 1.e-2;

% Finally, here are our SOS constraints, as defined above
F = sos(V-eps*(x'*x));
F = F + [sos(s)];
F = F + [sos(-nablaV*f-s*g-eps*(x'*x))];

% Now run the magic!
[sol, u, Q] = solvesos(F,[],[],[c;cs]);

if ~sol.problem
    disp("SOS Solution Found!")
    Vn = clean(u{1}'*Q{1}*u{1},1.e-6);
    sdisplay(Vn);
else
    disp("Constraints Failed. :(")
end
  
