function [t, x_traj] = rk42(x_dot_fun, tspan, x0)
% Runge Kutta method with fixed step size, as described on Wikipedia page.
% Currently enforces a zoh on control input.
%   Inputs: x_dot_fun - dynamics function pointer
%           tspan - time vector
%           x0
    [m, n] = size(x0);
    % Make x0 a column vector
    if n > m, x0=x0'; end;
    h = tspan(2) - tspan(1);
    % Preallocate state trajectory
    num_steps = length(tspan);
    x_traj = [x0'; zeros(num_steps-1, length(x0))];
    for n = 1:num_steps
%         disp(['t: ' num2str(t(n))]);
        tn = tspan(n);
        xn = x_traj(n, :)';
        % Zero-order hold on control input
        k1 = x_dot_fun(tn, xn);
        k2 = x_dot_fun(tn + h/2, xn + k1*h/2);
        k3 = x_dot_fun(tn + h/2, xn + k2*h/2);
        k4 = x_dot_fun(tn + h, xn + h*k3);
        x_traj(n+1, :) = (xn + (h/6)*(k1 + 2*k2 + 2*k3 + k4))';
    end
    t = tspan;
end