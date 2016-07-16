function [t, x_traj, u_traj] = rk4(x_dot_fun, u_fun, tspan, x0, num_steps)
% Runge Kutta method with fixed step size, as described on Wikipedia page.
% Currently enforces a zoh on control input.
%   Inputs: xdot - dynamics function pointer
%           tspan - [t0 tf]
%           x0
%           num_steps - total number of points in integration
%           u(t, x) - input function
    [m, n] = size(x0);
    % Make x0 a column vector
    if n > m, x0=x0'; end;
    h = (tspan(2) - tspan(1))/num_steps;
    % Preallocate state trajectory
    x_traj = [x0'; zeros(num_steps-1, length(x0))];
    t(1) = tspan(1);
    for n = 1:num_steps
%         disp(['t: ' num2str(t(n))]);
        xn = x_traj(n, :)';

        % Zero-order hold on control input
        un = u_fun(t(n), xn);
        u_traj(n) = un;
        k1 = x_dot_fun(t(n), xn, un);
        k2 = x_dot_fun(t(n) + h/2, xn + k1*h/2, un);
        k3 = x_dot_fun(t(n) + h/2, xn + k2*h/2, un);
        k4 = x_dot_fun(t(n) + h, xn + h*k3, un);

%         u_traj(n) = u_fun(t(n), xn);
%         k1 = x_dot_fun(t(n), xn, u_fun(t(n), xn) );
%         k2 = x_dot_fun(t(n) + h/2, xn + k1*h/2, u_fun(t(n)+h/2, xn) );
%         k3 = x_dot_fun(t(n) + h/2, xn + k2*h/2, u_fun(t(n)+h/2, xn) );
%         k4 = x_dot_fun(t(n) + h, xn + h*k3, u_fun(t(n)+h, xn) );
        
        x_traj(n+1, :) = (xn + (h/6)*(k1 + 2*k2 + 2*k3 + k4))';
        t(n+1) = t(n) + h;
    end
    t=t';
end