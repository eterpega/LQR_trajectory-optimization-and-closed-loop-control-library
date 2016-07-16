function [lqr, u_cl_fun, tIdxFun] = tvLqr(sys, lqr, tspan, x0, u0)
    % sys, sys.x_dot_sym, sys.x_dot_fun, sys.stateVars=[q1, q2, ... q1_dot, q2_dot, ...], sys.nStates
    % lqr.Qf, lqr.Q, lqr.S, lqr.nSteps
    
    t0 = tspan(1);
    tf = tspan(2);
    % Create Chebyshev polynomial representation of trajectories
    x0_p = chebfun(x0', 'equi');
    u0_p = chebfun(u0', 'equi');
    % Anonymous function for mapping time to range -1:1 for chebfun
    % representations
    tIdxFun = @(t) -1+2*(t-t0)/(tf-t0);
    % Create lists of system physical property names and values
    paramNames = fieldnames(sys.param);
    for n=1:length(paramNames), paramVals(n) = getfield(sys.param, paramNames{n}); end %#ok<AGROW,GFLD>    
    % Create symbolic x_dot function substituting in physical parameters
    physSys = subs(sys.x_dot_sym, paramNames, paramVals');
    % Linearize system
    A_lin_sym = jacobian(physSys, sys.stateVars); % df/dx
    B_lin_sym = jacobian(physSys, sys.inputVars(1));  % df/du
    
    A_lin_fun = matlabFunction(A_lin_sym, 'Vars', {[sys.stateVars; sys.inputVars]});
    B_lin_fun = matlabFunction(B_lin_sym, 'Vars', {[sys.stateVars; sys.inputVars]});
    
    % Create time-dependent (numeric) expressions for linearized A & B matrices
    % - these are actually dependent on the nominal trajectory state
    Alin_t = @(t) A_lin_fun([x0_p(tIdxFun(t))'; u0_p(tIdxFun(t))]);
    Blin_t = @(t) B_lin_fun([x0_p(tIdxFun(t))'; u0_p(tIdxFun(t))]);

    R_inv = inv(lqr.R);
    S_f = lqr.Q_f;
    S_dot = @(t, S, u) -(S*Alin_t(t) + Alin_t(t)'*S - S*Blin_t(t)*R_inv*Blin_t(t)'*S + lqr.Q);
    S_dot_wrapper = @(t, S, u) reshape(S_dot(t, reshape(S, sys.nStates, sys.nStates), u), sys.nStates^2, 1);
    % Solve differential Ricatti equation
    [S_t_traj, S_traj, ~] = rk4(S_dot_wrapper, @(t, S) 0, [tf t0], reshape(S_f, sys.nStates^2, 1), lqr.nSteps);
    K = zeros(length(S_t_traj), sys.nStates);
    eVals = zeros(length(S_t_traj), sys.nStates);
    for n=1:length(S_t_traj)
        Sn = reshape(S_traj(n, :), sys.nStates, sys.nStates);
        tn = S_t_traj(n);
        K(n, :) = R_inv*Blin_t(tn)'*Sn;
        eVals(n, :) = eig(Alin_t(tn)-Blin_t(tn)*K(n, :));
    end
    % Flip K so it's forwards in time
    lqr.K = flip(K, 1);
    lqr.eVals = flip(eVals, 1);
    lqr.K_p = chebfun(lqr.K, 'equi');
    u_cl_fun = @u_lqr;
%     u_ol_fun = @u_ol;
end

% V2 
function u = u_lqr(t, x, x0, u0, K, tIdx)
    % Expects K to be a chebfun polynomial representation of lqr gains
    % Same for x0 (nominal state trajectory) and u0 (nominal force trajectory) 
    u_nom = u0(tIdx(t));
    x_nom = x0(tIdx(t));
    Kn = K(tIdx(t));
    x_bar = x - x_nom';
    u_bar = -Kn*x_bar;
    u = u_bar + u_nom;
%     disp(['K: ' num2str(Kn) ', u0: ' num2str(u_nom), ', u~: ' num2str(u_bar) ', x~: ' num2str(x_bar')]);
%     fprintf('\nK: %f, u0: %f, u~: %f, x~: %f', Kn, u_nom, u_tilde, x_tilde);
end

% V1 - gain vector for each step of lqr simulation.
% function u = u_lqr(t, x, h, x0, u0, K, tIdx)
%     % Currently this controller requires that the rk4 simulation of the
%     % closed loop system be performed with the same size time steps as
%     % tvLqr(), i.e. a K gain vector is created for each step of the later
%     % CL simulation.
%     nStates = length(x0);
%     u_nom = u0(tIdx(t));
%     x_nom = x0(tIdx(t));
% %     x_nom = [];
% %     for k = 1:nStates, x_nom = [x_nom; x0{k}(tIdx(t))]; end
%     Kn = K(round(t/h)+1, :);
%     x_bar = x - x_nom';
%     u_bar = -Kn*x_bar;
%     u = u_bar + u_nom;
% %     disp(['K: ' num2str(Kn) ', u0: ' num2str(u_nom), ', u~: ' num2str(u_bar) ', x~: ' num2str(x_bar')]);
% %     fprintf('\nK: %f, u0: %f, u~: %f, x~: %f', Kn, u_nom, u_tilde, x_tilde);
% end

% function u = u_ol(t, x, h, u0)
%     if t == 0 
%         u = u0(1);
%     elseif rem(t, h) < .00001
%         % We're at a time which is a whole multiple of h
%         u = u0(round(t/h) + 1);
%     elseif rem(t*2, h) < .00001
%         % We're at a midpoint (knot point)
%         n = round((t-h)/h + 1);
%         % First order hold
%         u = (u0(n)+u0(n+1))/2;
%     else
%         error('Not a whole or half multiple of h');
%     end
% end

% function u = u_ol(t, x, h, u0)
%     u = u0(floor(t/h) + 1);
% end