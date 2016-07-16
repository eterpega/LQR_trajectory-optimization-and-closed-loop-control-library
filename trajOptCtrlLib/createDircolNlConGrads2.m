function sys = createDircolNlConGrads2(sys)
    % Creates analytic direct collocation gradients for a dynamical system.
    % Outputs a pointer to a function 
    %   gradf(x0, x1, u0, u1, T, t0, np)

    syms x0 x1 u0 u1 t0 T h np dummy real
    d = sym('d_%d', [1 sys.nStates]).';
    x0 = sym('x0_%d', [1 sys.nStates]).';
    x1 = sym('x1_%d', [1 sys.nStates]).';

    % Grab the eom for the system
    dynFun = sys.x_dot_sym.';
    % Create a list of variables that represent physical parameters of
    % system
    vars = symvar(dynFun);
    physPropVars = [];  % Physical properties variables (mass, length etc)
    for k = 1:length(vars)
        if isfield(sys.param, char(vars(k)))
            physPropVars = [physPropVars vars(k)];
        end
    end
    physPropVars = sort(physPropVars);
    % Swap qi, qi_dot etc for array of state values
    % Assumption is all coordinate variables begin with 'q' - need to fix
    % this
    for k = 1:sys.nStates/2
        dynFun = subs(dynFun, ['q' num2str(k) '_dot'], d(k+sys.nStates/2));
    end
    for k = 1:sys.nStates/2
        dynFun = subs(dynFun, ['q' num2str(k)], d(k));
    end
    % Create collocation expressions
    disp('Creating symbolic collocation expressions');
    f0 = subs(dynFun, ['t' 'u' d], [t0 u0 x0]);
    f1 = subs(dynFun, ['t' 'u' d], [t0+h u1 x1]);
    x_c = .5*(x0 + x1) + .125*h*(f0-f1);
    x_dot_c = 3*(x1 - x0)/(2*h) - .25*(f0+f1); 
    u_c = (u0 + u1)/2;  % First order hold (linear)
    t_c = (t0 + t0+h)/2;
    f_c = subs(dynFun, ['t' d 'u'], [t_c x_c u_c]);
    % Create defect expression
    defect = x_dot_c - f_c;
    defect = subs(defect, h, T/(np-1));
    assume(defect, 'real'); % Not sure if this helps, but can't hurt
    decvars = [x0; x1; u0; u1; T];
    disp('Calculating symbolic gradient expression');
    gradfSym = jacobian(defect, decvars).';
    disp('Saving numeric matlabFunction');
    % Save numeric function. One extra dummy argument so calling format
    % matches numeric central difference function
    sys.dcGrads = matlabFunction(gradfSym, 'File', [sys.name 'DcGrads'], 'vars', {x0, x1, u0, u1, T, t0, np, physPropVars, dummy});
end