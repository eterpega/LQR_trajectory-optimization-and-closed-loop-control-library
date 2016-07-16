    T = decVars(param.tIdx);
    x = reshape(decVars(param.xIdx(1):param.xIdx(2)), param.nStates, param.nKnotPoints);
    u = decVars(param.uIdx(1):param.uIdx(2));
    % Plot phase diagram
%     plotPhaseDiagram(x, u, param);
    % Create time vector
    t = linspace(0, T, param.nKnotPoints);
    % Step size
    h = t(2) - t(1);
    defects = zeros(param.nStates, param.nKnotPoints-1);
    for n = 1:param.nKnotPoints-1;
        x0 = x(:, n);
        x1 = x(:, n+1);
        u0 = u(n);
        u1 = u(n+1);
        f0 = param.dynFun(t(n), x0, u0, param);
        f1 = param.dynFun(t(n+1), x1, u1, param); % pendulum_cart_eom
        x_c = .5*(x0 + x1) + .125*h*(f0-f1);
        x_dot_c = 3*(x1 - x0)/(2*h) - .25*(f0+f1); 
        u_c = (u0 + u1)/2;  % First order hold (linear)
        t_c = (t(n)+t(n+1))/2;
        f_c = param.dynFun(t_c, x_c, u_c, param); 
        defects(:, n) = x_dot_c - f_c;
    end
    defects = reshape(defects, 1, param.nStates*(param.nKnotPoints-1));
%     boundaryEnd = (x(:, end) - param.xf)';
    ceq = [defects];% boundaryEnd];
    % No inequality constraints
    cineq = [];