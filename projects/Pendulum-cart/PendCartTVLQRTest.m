% TVLQR test - package together all the functions required to derive EOM
% for pendulum cart system, generate trajectory and create TVLQR

clear all;

% Generate EOM

% sysName = 'pendulumCartPointMass';
% syms m1 m2 l q1(t) q2(t) g b u
% coordVars = {q1, q2};
% % sys.nCoords = length(sys.coordVars);
% Q = [u; 0];
% T = .5*m1*diff(q1, 't')^2 + .5*m2*(diff(q1, 't') + l*diff(q2, 't')*cos(q2))^2 + ...
%     .5*m2*(l*diff(q2, 't')*sin(q2))^2;
% U = m2*g*l*(1-cos(q2));
% L = T-U;
% D = .5*b*diff(q1, 't')^2;
% disp('Deriving eom');
% % [sys.x_dot_sym, sys.x_dot_fun, sys.stateVars] = deriveEom(sysName, sys.coordVars, L, D, Q);
% % sys.nStates = length(sys.stateVars);
% sys = deriveEom(sysName, coordVars, L, D, Q);
load('pendulumCartPointMassSys.mat', 'sys');
syms u;
sys.inputVars = u; % Needed for LQR. Need to resolve better.
% % Need to think about this:
% sys.x_dot_fun = str2func(sys.x_dot_fun);

% Define physical properties of system
param.g = 9.81;
param.l = .5;
param.m1 = .5;
param.m2 = .5;
param.b = 0;
sys.param = param;

disp('Loading nominal trajectory');
% pendulumCartDirCol_120_ip_usq_200tsq_rk4
% [xnom, unom, T, param, tmp] = loadTrajectory('pendulumCartDirCol_120_ip_usq_200tsq_rk4'); % Works!
[xnom, unom, T, param, tmp] = loadTrajectory('pendulumCartPointMass_80_dircol_1usq_100uMx_(2)');
% [xnom, unom, T, param, tmp] = loadTrajectory('pendulumCartPointMass_100_dircol_1Tsq_0_005usq_60uMx'); % Works

% sys.param = param
[~, nKnotPoints] = size(xnom);
h = T/(nKnotPoints-1);
% Nominal trajectory time vector
t0 = linspace(0, T, nKnotPoints);

% Create LQR structure
% Gains for 120 knot point trajectory
lqr.Q = diag([10 20 1 1]);
lqr.R = 1*1;
lqr.Q_f = diag([10 20 1 1]); 
% lqr.Q_f = diag([20 20 1 1]); % WORKS
% % Gains for 60-knot point trajectory
% lqr.Q = diag([1 20 1 1]);
% lqr.R = 1*1;
% lqr.Q_f = diag([1 20 1 1]);

lqr.nSteps = nKnotPoints; % Create a gain for every knot point
% Get time varying LQR controller
disp('Calculating finite horizon LQR gains');
% [lqr, u_cl_fun, x0_p, u0_p, tIdxFun] = tvLqr(sys, lqr, [0 T], xnom, unom);
[lqr, u_cl_fun, tIdxFun] = tvLqrDirCol(sys, lqr, [0 T], xnom, unom);

% Change initial state
x_zero = [-1.5 -60*pi/180 0 -1*pi]'; % WORKS
% Perturb system physical properties
sys.param.l = .75;
% sys.param.m1 = .5;
% sys.param.m2 = .5;
sys.param.b = 1;
% ZOH input function
u_ol = @(t, x, h) unom(round(t/h) + 1);
% u_ol = @(t, x) xuOft(
disp('Simulating open- and closed-loop trajectories of perturbed system');
% Open loop simulation with as many steps as nominal trajectory
[t_vect_ol, x_traj_ol, u_traj_ol] = rk4(@(t, x, u) sys.x_dot_fun(t, x, u, sys.param), @(t, x) u_ol(t, x, h), [t0(1) t0(end)], x_zero, nKnotPoints-1);
% Closed loop simulation
% [t_vect, x_traj_cl, u_traj_cl] = rk4(@(t, x, u) dynFun(t, x, u, sys.param), @(t, x) u_cl_fun(t, x, h, x0_p, u0_p, lqr.K, tIdxFun), [0 T], x_zero, lqr.nSteps);
% [t_vect_cl, x_traj_cl, u_traj_cl] = rk4(@(t, x, u) sys.x_dot_fun(t, x, u, sys.param), @(t, x) u_cl_fun(t, x, x0_p, u0_p, lqr.K_p, tIdxFun), [0 T], x_zero, lqr.nSteps*10);
[t_vect_cl, x_traj_cl, u_traj_cl] = rk4(@(t, x, u) sys.x_dot_fun(t, x, u, sys.param), u_cl_fun, [0 T], x_zero, lqr.nSteps*10);

% [~, x_traj_ol, u_traj_ol] = rk4(@(t, x, u) dynFun(t, x, u, sys.param), @(t, x) u_ol_fun(t, sim.h, u0), [0 sim.tf], x_zero, sim.nSteps);

disp('Plot system response comparison');
% Plot comparison of state trajectories
plotTrajComp({t0, t_vect_ol, t_vect_cl}, {xnom, x_traj_ol, x_traj_cl}, 2, 2, ...
    [1 2 3 4], {':k', 'b', 'm'}, 'Pendulum cart (point mass)', ...
    {'q1', 'q2', 'q1 dot', 'q2 dot'}, {'Nominal', 'Open loop', 'Closed loop'});
% Input force trajectories
plotTrajComp({t0, t_vect_ol, t_vect_cl}, {unom, u_traj_ol, u_traj_cl}, 1, 1, ...
    [1], {':k', 'b', 'm'}, 'Pendulum cart (point mass)', ...
    {'u'}, {'Nominal', 'Open loop', 'Closed loop'}); %#ok<NBRAK>