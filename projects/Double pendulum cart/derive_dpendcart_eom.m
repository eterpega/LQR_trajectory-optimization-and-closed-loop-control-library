% Double pendulum with cart
clear variables;
syms m1 b1 m2 I2 c2 l2 b2 m3 I3 c3 l3 b3 q1(t) q2(t) q3(t) g u
sysName = 'doublePendCart';
coordVars = {q1, q2, q3};

x2 = q1+c2*sin(q2); y2 = -c2*cos(q2);
x3 = q1 + l2*sin(q2) + c3*sin(q3); y3 = -l2*cos(q2) - c3*cos(q3);
T = .5*m1*diff(q1, 't')^2 + .5*m2*diff(x2, 't')^2 + .5*m2*diff(y2, 't')^2 + ...
    .5*I2*diff(q2, 't')^2 + .5*m3*diff(x3, 't')^2 + .5*m3*diff(y3, 't')^2 + .5*I3*diff(q3, 't')^2;
V = m2*g*c2*(1 - cos(q2)) + m2*g*(l2*(1-cos(q2)) + c3*(1-cos(q3)));
L = T - V;
D = .5*b1*diff(q1, 't')^2 + .5*b2*diff(q2, 't')^2 + .5*b3*diff(q3, 't')^2;
Q = [u; 0; 0];
sys = deriveEom(sysName, coordVars, L, D, Q, true)
