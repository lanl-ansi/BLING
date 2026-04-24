function [x, par, history] = newton_solve(x, par, param, net, history, initial_solve)
 
x_old = x;
x = x(par.index.newton);

tol   = 10^(param.newton_error);
err   = 1e6;
new_it = 0;
delta = 0;

% Newton parameters
[C, par] = newton_parameters(par);

% Newton Solve
[x, delta] = solve(x, C, par, delta);

% Output data
output = [new_it; err];
x_old(par.index.newton) = x;
x = x_old;

% Incorporate Newton data in history
history = store_history(history, delta, output);

% Newton iteration with backtracking on residual
function [x, delta] = solve(x, C, par, delta)
    alpha_min = 1e-8;   % smallest allowed step length
    shrink    = 0.5;    % step reduction factor

    while err > tol
        % Residual and Jacobian at current iterate
        J = par.A_newton(x, 2, 3, 4) + C;
        F = (par.A_newton(x, 1, 1, 1) + C) * x - par.a_newton;

        % Newton step
        warning('off','all');
        dx = -(J \ F);
        warning('on','all');

        % Current residual norm
        err_old = max(abs(F));

        % Backtracking line search
        alpha = 1.0;
        while true
            x_trial = x + alpha * dx;
            F_trial = (par.A_newton(x_trial, 1, 1, 1) + C) * x_trial - par.a_newton;
            err_trial = max(abs(F_trial));

            % accept if residual decreases
            if err_trial < err_old
                break
            end

            alpha = alpha * shrink;

            if alpha < alpha_min
                break
            end
        end

        % Update iterate
        x     = x + alpha * dx;
        delta = max(abs(alpha * dx));
        err   = err_trial;
        new_it = new_it + 1;

        if new_it > param.newton_iter
            break
        end
    end
end

% --- Newton parameters ---
    function [C, par] = newton_parameters(par)
        % Compression ratio
        if initial_solve
            comp_ratio = net.comps(:,5);
        else
            comp_ratio = (par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes));
        end

        fr_node = net.comps(:,2);
        to_node = net.comps(:,3);

        C = spalloc(par.n.newton + par.n.slack, par.n.newton, 2*par.n.comps);
        row = par.n.slack + par.n.pipes + par.n.nodes + (1:par.n.comps);

        % Weighted incidence matrix
        C(row, par.index.nodes) = sparse([fr_node; to_node], ...
            [net.comps(:,1); net.comps(:,1)], ...
            [-comp_ratio; ones(par.n.comps,1)], ...
            par.n.nodes, par.n.comps)';

        % Newton Jacobian RHS
        par.a_newton = [par.a; zeros(par.n.comps,1)];
        row = par.n.slack + par.n.pipes + (1:par.n.nodes);
        if initial_solve
            par.a_newton(row) = par.q;
        else
            par.a_newton(row) = par.q ...
                - par.G * x_old(par.index.gnode_s) ...
                + par.G * x_old(par.index.gnode_d);
        end
    end

% --- Data storage ---
    function history = store_history(history, delta, output)
        k = numel(history) + 1;
        history(k).F_surp_val           = history(k-1).F_surp_val;
        history(k).F_comp_val           = history(k-1).F_comp_val;
        history(k).duality_gap          = history(k-1).duality_gap;
        history(k).delta                = delta;
        history(k).step                 = delta;
        history(k).constraint_violation = output(2);
        history(k).quality              = history(k-1).quality;
        history(k).lambda_gnode         = history(k-1).lambda_gnode;
        history(k).lambda_all           = history(k-1).lambda_all;
    end
end
