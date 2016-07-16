% Bunch of tests of the deriveEom function. 

clear variables;


% % Double pendulum
% syms L T U q1(t) q2(t) c1 c2 l1 l2 m1 I1 m2 I2 g t b1 b2 
% sysName = 'doublePendulum';
% coordVars = {q1, q2};
% % Pendulum 1
% xc1 = c1*sin(q1)
% yc1 = -c1*cos(q1)
% Ek1 = .5*m1*diff(xc1, 't')^2 + .5*m2*diff(yc1, 't')^2 + .5*I1*diff(q1, 't')^2
% Ep1 = m1*g*c1*(1-cos(q1))
% 
% % Pendulum 2
% xc2 = l1*sin(q1) + c2*sin(q2)
% yc2 = -l1*cos(q1) - c2*cos(q2)
% Ek2 = .5*m2*diff(xc2, 't')^2 + .5*m2*diff(yc2, 't')^2 + .5*I2*diff(q2, 't')^2
% Ep2 = m2*g*(l1*(1-cos(q1)) + c2*(1 - cos(q2)))
% 
% % Make Lagrangian
% L = Ek1 + Ek2 - (Ep1 + Ep2);
% D = .5*b1*diff(q1, 't')^2 + .5*b2*diff(q2, 't')^2;

% % Double pendulum with cart
% syms m1 b1 m2 I2 c2 l2 b2 m3 I3 c3 l3 b3 q1(t) q2(t) q3(t) g u
% sysName = 'doublePendCart';
% coordVars = {q1, q2, q3};
% 
% x2 = q1+c2*sin(q2); y2 = -c2*cos(q2);
% x3 = q1 + l2*sin(q2) + c3*sin(q3); y3 = -l2*cos(q2) - c3*cos(q3);
% T = .5*m1*diff(q1, 't')^2 + .5*m2*diff(x2, 't')^2 + .5*m2*diff(y2, 't')^2 + ...
%     .5*I2*diff(q2, 't')^2 + .5*m3*diff(x3, 't')^2 + .5*m3*diff(y3, 't')^2 + .5*I3*diff(q3, 't')^2;
% V = m2*g*c2*(1 - cos(q2)) + m2*g*(l2*(1-cos(q2)) + c3*(1-cos(q3)));
% L = T - V;
% D = .5*b1*diff(q1, 't')^2 + .5*b2*diff(q2, 't')^2 + .5*b3*diff(q3, 't')^2;
% Q = [u; 0; 0];
% sys = deriveEom(sysName, coordVars, L, D, Q)

% % Two masses with vertical wall, springs between masses and to wall
% syms m1 m2 k1 k2 q1(t) q2(t) b1 b2 u
% coordVars = {q1, q2};
% sysName = 'twoMassesSpringsWall';
% T = .5*m1*diff(q1, 't')^2 + .5*m2*diff(q2, 't')^2;
% U = + .5*k1*q1^2 + .5*k2*(q2 - q1)^2;
% D = 0;
% Q = [u; 0];
% L = T-U;
% [x_dot_sym, eomFile, stateVars] = deriveEom(sysName, coordVars, L, D, Q);


% % Simple pendulum
syms m1 l1 q1(t) g u real
sysName = 'pointMassPendulum';
coordVars = {q1};
T = .5*(m1*l1^2)*diff(q1, 't')^2;
V = m1*g*l1*(1-cos(q1));
L = T-V;
D = 0;
Q = u;
sys = deriveEom(sysName, coordVars, L, D, Q, true);

% % Two masses, springs on leftmost mass and between them
% sysName = 'doubleMassSprings';
% coordVars = {q1, q2};
% syms k1 k2
% T = .5*m1*diff(q1, 't')^2 + .5*m2*diff(q2, 't')^2;
% V = .5*k1*q1^2 + .5*k2*q2^2;
% L = T - V;
% D = 0;

% % Pendulum cart
% sysName = 'pendulumCartTest2';
% syms m1 m2 c1 q1(t) q2(t) g b u
% coordVars = {q1, q2};
% Q = [u; 0];
% T = .5*m1*diff(q1, 't')^2 + .5*m2*(q1 + c1*sin(diff(q2, 't')))^2 + .5*m2*(c1*cos(diff(q2, 't')))^2;
% U = m2*g*c1*(1-cos(q2));
% L = T-U;
% D = .5*b*diff(q1, 't')^2;
% %     D = 0; % No damping
% [q_ddot_symvar, q_ddot_funcname] = deriveEom(sysName, coordVars, L, D, Q);

% sysName = 'pendulumCartPointMass';
% syms t m1 m2 l q1(t) q2(t) g b u real
% sys.coordVars = {q1, q2};
% sys.nCoords = length(sys.coordVars);
% sys.inputVars = u; % Needed for LQR. Need to resolve better.
% Q = [u; 0];
% % T = .5*m1*diff(q1, 't')^2 + .5*m2*(q1 + l*sin(diff(q2, 't')))^2 + .5*m2*(l*cos(diff(q2, 't')))^2;
% T = .5*m1*diff(q1, 't')^2 + .5*m2*(diff(q1, 't') + l*diff(q2, 't')*cos(q2))^2 + ...
%     .5*m2*(l*diff(q2, 't')*sin(q2))^2;
% U = m2*g*l*(1-cos(q2));
% L = T-U;
% % Hand derivation
% % L = .5*(m1+m2)*diff(q1, t)^2 + .5*m2*l^2*diff(q2, t)^2 + diff(q1, t)*diff(q2, t)*m2*l*cos(q2) - m2*g*l*(1-cos(q2));
% 
% % D = .5*b*diff(q1, 't')^2;
% D = 0;
% sys = deriveEom(sysName, sys.coordVars, L, D, Q, true)

% % Book example
% sysName = 'pendulumCartTest';
% coordVars = {q1, q2};
% L = .5*m1*diff(q1, 't')^2 + .5*m2*((diff(q1, 't') - l2*diff(q2, 't')*sin(q2))^2 + (l2*diff(q2, 't')*cos(q2))^2) - m2*g*l2*sin(q2);
% D = 0;

% Wikipedia double pendulum example
% L = .5*m_1*l_1^2*(diff(q_2, 't')^2+4*diff(q_1, 't')^2 + 3*diff(q_1, 't')*diff(q_2, 't')*cos(q_1-q_2)) + .5*m_1*g*l_1*(3*cos(q_1)+cos(q_2));
