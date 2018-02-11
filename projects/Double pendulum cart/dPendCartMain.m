% Generate an optimised trajactory
clear variables;
% derive_dpendcart();
load('doublePendCartSys.mat', 'sys');

% % % % Parameters from SW model, 29/10/17, 400mm pendulums
% % % % Cart
% % % p.m1 = 0.244; p.b1 = 0;%0.01;
% % % % Pendulum 1
% % % p.m2 = 0.114; p.I2 = 0.011831921; p.c2 = 0.277; p.l2 = 0.432; p.b2 = 0.00125;
% % % % Pendulum 2
% % % p.m3 = 0.04253; p.I3 = 0.0030585; p.c3 = 0.1976; p.l3 = 0.420; p.b3 = 0.00125;
% % % p.g = 9.81;

% Parameters from SW model, 13/11/17, 300mm pendulums
% Cart
p.m1 = 0.244; p.b1 = 0;%0.01;
% Pendulum 1
p.m2 = 0.1017; p.I2 = 0.006585751; p.c2 = 0.2163; p.l2 = 0.3268; p.b2 = 0.00125;
% Pendulum 2
p.m3 = 0.0397; p.I3 = 0.001652098; p.c3 = 0.1499; p.l3 = 0.320; p.b3 = 0.00125;
p.g = 9.81;


sys.param = p;

% % If parameters change, run these two lines to update gradient function
% sys = createDircolNlConGrads2(sys);
% save([sys.name 'Sys.mat'], 'sys');
% return;

method = 'dircol';
gradients = 'centraldiff';  % solvergrads/centraldiff/analytic
 
nPoints = 250;
x0 = [0 0 0 0 0 0]';
xf = [0 pi pi 0 0 0]';
% x0 = [0 pi pi 0 0 pi/2]';
% xf = [0 pi pi 0 0 pi/2]';

% guess = '1005'; %'doublePendCart_80_dircol_1usq_50uMx';%'doublePendCart_30_dircol_1Tsq_1usq_40uMx';

% guess.traj = (xf-x0)*linspace(0, 1, nPoints);
% guess.T = 3;
% guess.traj(5:6,:) = (pi/guess.T)*ones(2, nPoints);
% guess.u = zeros(1, nPoints);

guess = 0;

xLims = [-0.6 -Inf -Inf -20 -Inf -Inf; 0.6 Inf Inf 20 Inf Inf]';
% xLims = [-0.6 -Inf -Inf -20 -Inf pi/2*.9; 0.6 Inf Inf 20 Inf pi/2*1.1]';
uMax = 50;
tLims = [0 5];
cost.u = 0;%.5;%30/nPoints;
cost.uSmooth = 0;%0.1;
cost.accSmooth = 0;%10;
cost.T = 1;%10;%50;
[traj, u, T, param, exitflag, output, duration] = trajOpt(sys, method, gradients, cost, nPoints, x0, xf, guess, xLims, uMax, tLims);