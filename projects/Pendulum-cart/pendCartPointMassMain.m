clear all;
load('pendulumCartPointMassSys.mat', 'sys');

% p.g = 9.81;
% p.l = .5;
% p.m1 = .5;
% p.m2 = .5;
% p.b = 0;
% sys.param = p;

% % If parameters change, run these two lines to update gradient function
% sys = createDircolNlConGrads2(sys);
% save([sys.name 'Sys.mat'], 'sys');
% return

p.m1 = 2.0; p.m2 = 0.5; p.g = 9.81; p.l = 0.5; % p.b = 0;
sys.param = p;

method = 'dircol';
gradients = 'analytic'; % solvergrads/centraldiff/analytic
nPoints = 80;
x0 = [0 0 0 0]';
% xf = [0.8 pi 0 0]';
xf = [0 pi 0 0]';

% guess.traj = (xf-x0)*linspace(0, 1, nPoints);
% guess.T = 3;
% guess.traj(4,:) = (pi/guess.T)*ones(1, nPoints);
% guess.u = zeros(1, nPoints);
guess = 0;
% guess = 'pendulumCartPointMass_40_dircol_1Tsq_20uMx_(2)';


xLims = [-2*.8 -2*pi -Inf -Inf; 2*.8 2*pi Inf Inf]';
% xLims = [-2*.8 -Inf -Inf -Inrf; 2*.8 Inf Inf Inf]';
uMax = 20;
tLims = [.1 5];
% tLims = [.1 5];
cost.u = 0;%1/nPoints;
cost.T = 1;
[traj, u, T, param, exitflag, output] = trajOpt2(sys, method, gradients, cost, nPoints, x0, xf, guess, xLims, uMax, tLims);

% % Construct name and save trajectory
% str = [sys.name '_' num2str(nPoints) '_' method '_'];
% if cost.T > 0, str = [str num2str(cost.T, 2) 'Tsq_']; end
% if cost.u > 0, str = [str num2str(cost.u, 2) 'usq_']; end
% str = [str num2str(uMax) 'uMx'];
% filename = str;
% condition = true;
% n = 2;
% while condition
%     if exist([filename '.mat'], 'file') == 2
%         filename = [str '_(' num2str(n) ')'];
%     else
%         condition = false;
%     end
%     n = n + 1;
% end
% filename = strrep(filename, '.', '_');
% disp(['Saving trajectory ' filename]);
% saveTrajectory(filename);