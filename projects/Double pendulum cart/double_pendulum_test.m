% Double pendulum test
load('doublePendCartSys.mat', 'sys');
% 10 swg thick ali tube
p.m1 = .125; p.b1 = .05;
p.m2 = 0.19931639; p.I2 = 0.00338952; p.c2 = 0.21172; p.l2 = 0.422; p.b2 = 0;%0.00125;
% 2nd pendulum w/o bob
p.m3 = 0.186976; p.I3 = 0.002831; p.c3 = 0.19704; p.l3 = 0.411; p.b3 = 0;%0.00125;
p.g = 9.81;

% p.m1 = 0.125; p.b1 = 0;%0.01;
% p.m2 = 0.098; p.I2 = 0.00205; p.c2 = 0.20956; p.l2 = 0.422; p.b2 = 0.00125;
% % 2nd pendulum with bob at end
% % p.m3 = 0.13971; p.I3 = 0.0030865; p.c3 = 0.26557; p.l3 = 0.4; p.b3 = 0.005;
% % 2nd pendulum w/o bob
% p.m3 = 0.08646; p.I3 = 0.001446; p.c3 = 0.18; p.l3 = 0.411; p.b3 = 0.00125;
% p.g = 9.81;

sys.param = p;

x0 = [0 -135*pi/180 0 0 0 0]';
[t, traj] = ode45(@(t, x) sys.x_dot_fun(t, x, 0, sys.param), [0 20], x0);
T = t(end);
traj = traj';
animDoublePendCart;