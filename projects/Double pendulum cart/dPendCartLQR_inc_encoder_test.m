% dPendCartLQR
clear variables
disp('Loading system model');
load('doublePendCartSys.mat', 'sys');
syms u;
sys.inputVars = u; % Needed for LQR. Need to resolve better.
% % 3/4 in x 18 swg
% p.m1 = 0.125; p.b1 = 0;%0.01;
% p.m2 = 0.098; p.I2 = 0.00205; p.c2 = 0.20956; p.l2 = 0.422; p.b2 = 0.00125;
% % 2nd pendulum with bob at end
% p.m3 = 0.13971; p.I3 = 0.0030865; p.c3 = 0.26557; p.l3 = 0.4; p.b3 = 0.005;
% % % 2nd pendulum w/o bob
% % p.m3 = 0.08646; p.I3 = 0.001446; p.c3 = 0.18; p.l3 = 0.411; p.b3 = 0.00125;
% p.g = 9.81;
% sys.param = p;

disp('Loading nominal trajectory');
% Max ang vel 12 rad/s: doublePendCart_150_dircol_1usq_50uMx
% Max ang vel 12 rad/s doublePendCart_300_dircol_1usq_50uMx <- one I've worked on a lot
% Max ang vel 16 rad/s: doublePendCart_150_dircol_1Tsq_1usq_50uMx
% 1st equilib. Max ang vel ~8 rad/s: doublePendCart_150_dircol_1usq_25uMx
% Smoothed u: doublePendCart_80_dircol_1usq_50uMx

% lqr.Q = .1*diag([1 5 5 0 0 0]); lqr.R = 1; lqr.Q_f = lqr.Q;
<<<<<<< HEAD:projects/Double pendulum cart/dPendCartLQR_inc_encoder_test.m
% doublePendCart_150_dircol_1usq_25uMx % First equilibrium point
[xnom, unom, T, param, tmp] = loadTrajectory('doublePendCart_150_dircol_1usq_50uMx');   %, doublePendCart_240_dircol_1Tsq_1usq_50uMx, doublePendCart_150_dircol_1usq_25uMx  doublePendCart_240_dircol_1Tsq_1usq_50uMx
=======

[xnom, unom, T, param, tmp] = loadTrajectory('1013');%doublePendCart_150_dircol_1usq_25uMx');   %, doublePendCart_240_dircol_1Tsq_1usq_50uMx, doublePendCart_150_dircol_1usq_25uMx  doublePendCart_240_dircol_1Tsq_1usq_50uMx
>>>>>>> 517516febf6516fe4b8dab865e58beea24697ad4:projects/Double pendulum cart/dPendCartLQR.m
sys.param = param.physProp;
%'doublePendCart_120_dircol_10Tsq_0_25usq_40uMx'); % Works
[~, nPoints] = size(xnom);
h = T/(nPoints-1);
% Nominal trajectory time vector
t0 = linspace(0, T, nPoints);

% Create LQR structure
lqr.Q = 1*eye(6);
% lqr.Q = .0001*diag([1 5 5 1 5 5]);
<<<<<<< HEAD:projects/Double pendulum cart/dPendCartLQR_inc_encoder_test.m
% lqr.Q = diag([1 2.5 2.5 .05 .1 .1]); % Works with doublePendCart_150_dircol_1usq_25uMx
=======
lqr.Q = diag([.5 2.5 2.5 .05 .1 .1]);
>>>>>>> 517516febf6516fe4b8dab865e58beea24697ad4:projects/Double pendulum cart/dPendCartLQR.m
lqr.R = 1;
% lqr.Q_f = .5*eye(6);%diag([5 5 5 1 1 1]); %5*eye(sys.nStates);
% lqr.Q_f = zeros(6, 6);%diag([1 5 5 .5 .5 .5]);
lqr.Q_f = lqr.Q;

lqr.nSteps = nPoints; % Create a gain for every knot point
% Get time varying LQR controller
disp('Calculating finite horizon LQR gains');
% [lqr, u_cl_fun, x0_p, u0_p, tIdxFun] = tvLqr(sys, lqr, [0 T], xnom, unom);
[lqr, u_cl_fun, tIdxFun] = tvLqrDirCol(sys, lqr, [0 T], xnom, unom);
% return;
% Change initial state
<<<<<<< HEAD:projects/Double pendulum cart/dPendCartLQR_inc_encoder_test.m
x_zero = [0 0 0 0 0 0]';
% x_zero = [-.5 -30*pi/180 30*pi/180 0 0 0]';
% % Perturb system physical properties
sys.param.b1 = 0.25;
% sys.param.b2 = 0.00075;
sys.param.m2 = sys.param.m2*1.05;
=======
x_zero = [0 0 0 0 0 0];
% x_zero = [-.5 -30*pi/180 30*pi/180 0 0 0]';
% % Perturb system physical properties
sys.param.b1 = 0.1;
sys.param.b2 = sys.param.b2*1.05;
>>>>>>> 517516febf6516fe4b8dab865e58beea24697ad4:projects/Double pendulum cart/dPendCartLQR.m
sys.param.b3 = sys.param.b3*1.05;
sys.param.m2 = sys.param.m2*1.05;
sys.param.m3 = sys.param.m3*1.05;

% sys.param.b
% ZOH input function
u_ol = @(t, x) unom(round(t/h) + 1);
% u_ol = @(t, x) xuOft(
disp('Simulating open- and closed-loop trajectories of perturbed system');

% Open loop simulation with as many steps as nominal trajectory
% [t_vect_ol, x_traj_ol, u_traj_ol] = rk4(@(t, x, u) sys.x_dot_fun(t, x, u, sys.param), u_ol, [t0(1) t0(end)], x_zero, nKnotPoints-1);
ol_sol = ode45(@(t, x) sys.x_dot_fun(t, x, u_ol(t, x), sys.param), [t0(1) t0(end)], x_zero);
x_traj_ol = deval(ol_sol, t0);
 
% Initialise incremental encoder
pulleyRad = 0.03056/2;  % Cart pulley radius (m)
% incEnc(0, x_zero, 100, [1440*(1/pulleyRad)/(2*pi) 1440 1440]');
enc = IncEnc(3, 0, 100, (80000/1440)*[1440*(1/pulleyRad)/(2*pi) 1440 1440]', true);
% enc = IncEnc(3, 0, 100, 0, true);

% Closed loop simulation
% [t_vect_cl, x_traj_cl, u_traj_cl] = rk4(@(t, x, u) sys.x_dot_fun(t, x, u, sys.param), u_cl_fun, [0 T], x_zero, lqr.nSteps*10);
% cl_sol = ode45(@(t, x) sys.x_dot_fun(t, x, u_cl_fun(t, x), sys.param), [t0(1) t0(end)], x_zero);
% x_traj_cl = deval(cl_sol, t0);
% t_cl = t0;
% Fixed step simulation

[t_cl, x_traj_cl] = rk42(@(t, x) sys.x_dot_fun(t, x, u_cl_fun(t, enc.read(t, x)), sys.param), t0(1):0.001:t0(end), x_zero);
% [t_cl, x_traj_cl] = rk42(@(t, x) sys.x_dot_fun(t, x, u_cl_fun(t, x), sys.param), t0(1):0.001:t0(end), x_zero);

disp('Plot system response comparison');
% Plot comparison of state trajectories
plotTrajComp({t0, t0, t_cl}, {xnom, x_traj_ol, x_traj_cl}, 2, 3, ...
    [1 2 3 4 5 6], {':k', 'b', 'm'}, 'Double pendulum cart', ...
    {'q1', 'q2', 'q3', 'q1 dot', 'q2 dot', 'q3 dot'}, {'Nominal', 'Open loop', 'Closed loop'});
% plotTrajComp({t0, t0, t0}, {xnom, x_traj_ol, x_traj_cl}, 2, 3, ...
%     [1 2 3 4 5 6], {':k', 'b', 'm'}, 'Double pendulum cart', ...
%     {'q1', 'q2', 'q3', 'q1 dot', 'q2 dot', 'q3 dot'}, {'Nominal', 'Open loop', 'Closed loop'});


% Plot encoder estimate of state vs true state:
plotTrajComp({enc.Trajectory.xIn.t, enc.Trajectory.xEst.t}, {enc.Trajectory.xIn.traj, enc.Trajectory.xEst.traj}, 2, 3, ...
    [1 2 3 4 5 6], {'b', 'm--'}, 'Double pendulum cart', ...
    {'q1', 'q2', 'q3', 'q1 dot', 'q2 dot', 'q3 dot'}, {'Actual', 'Estimated'});

% Input force trajectories
% plotTrajComp({t0, t_vect_ol, t_vect_cl}, {unom, u_traj_ol, u_traj_cl}, 1, 1, ...
%     [1], {':k', 'b', 'm'}, 'Pendulum cart (point mass)', ...
%     {'u'}, {'Nominal', 'Open loop', 'Closed loop'}); %#ok<NBRAK>
<<<<<<< HEAD:projects/Double pendulum cart/dPendCartLQR_inc_encoder_test.m
% % plotTrajComp({t0, t0}, {unom, u_traj_cl}, 1, 1, ...
% %     [1], {'b', 'm'}, 'Pendulum cart (point mass)', ...
% %     {'u'}, {'Open loop', 'Closed loop'}); %#ok<NBRAK>
=======
plotTrajComp({t0, t0}, {unom, u_traj_cl}, 1, 1, ...
    [1], {'b', 'm'}, 'Pendulum cart (point mass)', ...
    {'u'}, {'Open loop', 'Closed loop'}); %#ok<NBRAK>
>>>>>>> 517516febf6516fe4b8dab865e58beea24697ad4:projects/Double pendulum cart/dPendCartLQR.m
