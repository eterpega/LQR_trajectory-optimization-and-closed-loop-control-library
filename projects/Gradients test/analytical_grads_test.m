% Analytical gradients test
% Sets up a trajectory optimization problem but just for testing 
clear variables

load('pointMassPendulumSys.mat')
p.m1 = 2; p.l1 = 1; p.g = 9.81;
sys.param = p;

% % If parameters change, run these two lines to update gradient function
% sys = createDircolNlConGrads2(sys);
% save([sys.name 'Sys.mat'], 'sys');

method = 'dircol';
nPoints = 20;
x0 = [0 0]';
xf = [pi 0]';

guess.traj = (xf-x0)*linspace(0, 1, nPoints);
guess.T = 1.5;
% guess.traj(4:5,:) = (pi/guess.T)*ones(2, nPoints);
guess.u = zeros(1, nPoints);

xLims = [-Inf -Inf; Inf Inf]';
uMax = 25;
tLims = [.1 3];
cost.u = 0;%30/nPoints;
cost.T = 1;%10;%50;
gradients = 'analytic'; % solvergrads/centraldiff/analytic
[traj, u, T, param, exitflag, output] = trajOpt2(sys, method, gradients, cost, nPoints, x0, xf, guess, xLims, uMax, tLims);
t = linspace(0,T, length(traj));
figure
plot(t, traj(1,:))