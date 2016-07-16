% TVLQR test

% Simple pendulum. Create trajectory from abitrary (sinusoidal) input
% function.
% Same as simpleLinTrajTvLQR(), except linearized A and B functions are
% represented using Chebyshev polynomial with Chebfun (simplifies things a bit).

syms th th_dot u t
% sym(t)

m = 1;
l = 1;
b = 1;
g = 9.81;
t0 = 0;
tf = 10;
xz = [0; 0];
nSteps = 100;
% sim.h = (sim.tf-sim.t0)/sim.nSteps;

sys.nStates = 2;
sys.stateVars = {th; th_dot};
sys.inputVars = [u];
% Pendulum eom (symbolic)
sys.x_dot_sym = [th_dot; u/(m*l^2) - (g/l)*sin(th) - b*th_dot/(m*l^2)];
% Numeric version
stateVars = [sys.stateVars{1}; sys.stateVars{2}];
sys.x_dot_fun = matlabFunction(sys.x_dot_sym, 'vars', {t, stateVars, sys.inputVars});

% Abitrary input function
u_f = @(t, x) sin(.2*2*pi*t);
% Simulate system and get trajectory
[t_vect, x0, u0] = rk4(sys.x_dot_fun, u_f, [t0 tf], xz, nSteps);

% Create LQR structure
lqr.Q = 5*eye(sys.nStates);
lqr.R = .2*1;
lqr.Q_f = 2*eye(2);
% Create sim structure
sim.t0 = t0;
sim.tf = tf;
sim.nSteps = 200;
sim.h = (sim.tf-sim.t0)/sim.nSteps;
% Get time varying LQR controller
[lqr, u_cl_fun, x0_p, u0_p, tIdxFun] = tvLqr(sys, sim, lqr, x0, u0);

% Calculate open loop of perturbed system w. nominal input
b = 2.5;
l = 1.05;
m = .95;
xz = [-35*pi/180; .5];

% Plot the nominal trajectory
figure;
subplot(2, 1, 1); hold on; grid on;
ylabel('theta (rad)');
plot(t_vect, x0(:, 1), '-.', 'LineWidth', 1, 'Color', [.5 .5 .5]);
subplot(2, 1, 2); hold on; grid on; ylabel('theta dot (rad/s)'); xlabel('t (s)'); 
plot(t_vect, x0(:, 2), '-.', 'LineWidth', 1, 'Color', [.5 .5 .5]);

q_dot_fun2 = matlabFunction(sys.x_dot_sym, 'vars', {t, stateVars, sys.inputVars});
u_ol_fun = @(t, ~) u0_p(tIdxFun(t));
[t_vect, x_traj_ol, ~] = rk4(q_dot_fun2, u_ol_fun, [t0 tf], ...
    xz, sim.nSteps);
subplot(2,1,1);
plot(t_vect, x_traj_ol(:, 1), 'Color', [.7 .9 .7], 'LineWidth', 1);
subplot(2, 1, 2);
plot(t_vect, x_traj_ol(:, 2), 'Color', [.7 .9 .7], 'LineWidth', 1);

% Closed loop trajectory of perturbed system
[t_vect, x_traj_cl, u_traj_cl] = rk4(q_dot_fun2, @(t, x) ...
    u_cl_fun(t, x, sim.h, x0_p, u0_p, lqr.K, tIdxFun), [t0 tf], ...
    xz, sim.nSteps);
subplot(2, 1, 1);
plot(t_vect, x_traj_cl(:, 1), 'Color', [.5 .5 1], 'LineWidth', 1);
subplot(2, 1, 2);
plot(t_vect, x_traj_cl(:, 2), 'Color', [.5 .5 1], 'LineWidth', 1);
legend('Nominal trajectory', 'Perturbed system, open loop', 'Perturbed system, closed loop');    


