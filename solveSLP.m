function [x, history, par] = solveSLP(net, par, param, x, history)
% Successive Linear Program with line search + trust region
%
% Key changes:
% 1) Trust-region update uses the accepted step s = alpha * dx_full
% 2) Trust-region ratio uses actual/predicted reduction:
%       ared = Psi_curr - Psi_trial
%       pred = Psi_curr - Psi_model
%       rho_tr = ared / pred
% 3) Armijo test is based on predicted reduction from the same model
% 4) Trust-region expansion only occurs when the accepted step is near the boundary

options = optimoptions('linprog','Display','none','Algorithm','dual-simplex-highs');

% Armijo backstepping defaults
if ~isfield(param, 'armijo_eta'),    param.armijo_eta    = 1e-6; end
if ~isfield(param, 'backstep_beta'), param.backstep_beta = 0.5;  end
if ~isfield(param, 'alpha_min'),     param.alpha_min     = 1e-2; end
if ~isfield(param, 'rho'),           param.rho           = 1e-3; end

eta       = param.armijo_eta;
beta_ls   = param.backstep_beta;
alpha_min = param.alpha_min;
rho       = param.rho;

LMP_gnode  = par.n.slack + par.n.pipes + par.G' * par.index.nodes;
LMP_all    = par.n.slack + par.n.pipes + par.index.nodes;
comp_price = net.comps(:,8);

delta = param.step_size;
x_new = x;

% Local compressor energy objective
[par, Fval] = SS_exact_comp(par, param);

% Precompute constants used in history metrics
den_surp = sum(par.xmax(par.index.gnode_s));
den_qual = sum(par.xmax(par.index.gnode_d) + par.G' * par.q);
Gtq      = par.G' * par.q;

% Solve optimization with SLP
[x, history] = slp(x, history);

% Remove NAN entries
cons_viol = [history.constraint_violation];
valid = ~isnan(cons_viol);
history = history(valid);

% --- Function definitions ---

    function [x, history] = slp(x, history)

        r_curr   = -par.a + par.A(x,1,1,1) * x;
        Phi_curr = NLP_obj(x, par);
        Psi_curr = merit(r_curr, Phi_curr, rho);

        for k = 1:param.max_iter

            % Linearized equality constraints: J * dx = -r_curr
            J = par.A(x,2,3,4);

            % Local compressor energy
            if param.exact_comp_energy == 1
                par = SS_exact_comp_calculation(x, par);
            end

            % Bounds on dx
            dxmin = max(-delta, par.xmin - x);
            dxmax = min( delta, par.xmax - x);

            % Local linear solve
            [dx_full, ~, exitflag, ~, lambda] = linprog( ...
                par.L', par.B, par.b - par.B * x, J, -r_curr, dxmin, dxmax, options);

            % LP infeasible
            if exitflag ~= 1
                dx_fail = zeros(size(x));
                [history, ~] = store_history(par, exitflag, k, x_new, dx_fail, history, delta);
                if k < 10
                    delta = delta + 0.5 * param.step_size;
                else
                    delta = delta + 0.1 * param.step_size;
                end
                continue
            end

            % Predicted reduction for full step, useful for diagnostics
            s_full         = dx_full;
            Psi_model_full = model_merit(r_curr, J, s_full, Phi_curr, par, rho);
            pred_full      = Psi_curr - Psi_model_full;

            % Armijo backstepping line search on the merit function
            alpha      = 1;
            accepted   = false;
            x_trial    = x;
            r_trial    = r_curr;
            Phi_trial  = Phi_curr;
            Psi_trial  = Psi_curr;
            Psi_model  = Psi_curr;
            pred       = 0;

            while alpha >= alpha_min
                s = alpha * dx_full;
                x_trial = x + s;

                r_trial   = -par.a + par.A(x_trial,1,1,1) * x_trial;
                Phi_trial = NLP_obj(x_trial, par);
                Psi_trial = merit(r_trial, Phi_trial, rho);

                % Merit model evaluated at accepted step s
                Psi_model = model_merit(r_curr, J, s, Phi_curr, par, rho);
                pred      = Psi_curr - Psi_model;

                % Require positive predicted reduction
                if pred <= 0
                    alpha = beta_ls * alpha;
                    continue
                end

                % Armijo condition based on model-predicted reduction
                if Psi_trial <= Psi_curr - eta * pred
                    accepted = true;
                    break
                end

                alpha = beta_ls * alpha;
            end

            if ~accepted
                % Step rejected by line search: shrink trust region
                delta = max(0.5 * delta, 1e-3);
                continue
            end

            % Accepted step
            dx    = alpha * dx_full;
            x_new = x_trial;
            r_new = r_trial;

            % Actual/predicted reduction ratio
            ared   = Psi_curr - Psi_trial;
            rho_tr = ared / (pred + eps);

            % Store history
            [history, ~] = store_history(par, exitflag, k, x_new, dx, history, delta, lambda, r_new);

            % Update trust region using accepted step and TR ratio
            delta = trust_region(param, dx, delta, rho_tr);

            % Update iterate
            x        = x_new;
            r_curr   = r_new;
            Phi_curr = Phi_trial;
            Psi_prev = Psi_curr;
            Psi_curr = Psi_trial;

            % Check convergence
            if abs(Psi_prev - Psi_curr) < 10^(param.nom_quality_tol)
                if max(abs(r_curr)) < 10^(param.error_tol)
                    break
                end
            end

            fprintf(['Iter: %d | Residual: %.3e | Merit: %.3e | Model Merit: %.3e | ' ...
                     'ared: %.3e | pred: %.3e | rho_tr: %.3e | alpha: %.3e | delta: %.3e\n'], ...
                     k, norm(r_new, inf), Psi_curr, Psi_model, ared, pred, rho_tr, alpha, delta);
        end
    end

% Data storage
    function [history, r_new] = store_history(par, exitflag, k, x_new, dx, history, delta, lambda, r_new)
        if exitflag ~= 1 && k == 1
            history(k).F_comp_val = 0.5;
            history(k).F_surp_val = 1;
            history(k).duality_gap = 0.5;
            history(k).delta = delta;
            history(k).step = delta;
            history(k).constraint_violation = 0.5;
            history(k).lambda_gnode = zeros(length(LMP_gnode),1);
            history(k).lambda_all = zeros(length(LMP_all),1);
            history(k).merit = 0.5;
            r_new = zeros(size(par.a));
            history(k).quality = 0;
            return
        elseif exitflag ~= 1 && k > 1
            history(k) = history(k-1);
            r_new = zeros(size(par.a));
            return
        end

        if exitflag == 1
            cons_violation = norm(r_new, inf);
            F_new = Fval(x_new);
            F_old = Fval(x_new - dx);
            duality_gap = norm( ...
                F_new - F_old ...
                - 1/(1-param.Econ_weight) * par.L(1:par.index.pipes(end)) * dx(1:par.index.pipes(end)), ...
                inf);

            history(k).F_surp_val = sum(x_new(par.index.gnode_s)) / den_surp;
            history(k).F_comp_val = F_new;
            history(k).duality_gap = max(duality_gap, 1e-12);
            history(k).delta = delta;
            history(k).step = max(abs(dx));
            history(k).constraint_violation = cons_violation;
            history(k).quality = sum(x_new(par.index.gnode_d) + Gtq) / den_qual;
            history(k).lambda_gnode = -lambda.eqlin(LMP_gnode);
            history(k).lambda_all   = -lambda.eqlin(LMP_all);
            history(k).merit = cons_violation + (1 - param.Econ_weight) * F_new;
        end
    end

    function [par, Fval] = SS_exact_comp(par, param)
        r = param.specific_heat_capacity_ratio;
        ratio = (r - 1) / r;

        if param.exact_comp_energy == 1
            a = (1 - param.Econ_weight) * (1 - r) / r;
            b = (1 - param.Econ_weight);
            par.alpha = @(x) a * comp_price .* x(par.index.comps)./(par.Cp * x(par.index.nodes)) .* ...
                ((par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes))) .^ ratio;
            par.beta  = @(x) a * comp_price .* x(par.index.comps)./(par.Cm * x(par.index.nodes)) .* ...
                ((par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes))) .^ ratio;
            par.gamma = @(x) b * max(comp_price .* ...
                (((par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes))) .^ ratio - 1), 0);
        end

        Fval = @(x) sum(comp_price .* x(par.index.comps) .* ...
            ((max((par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes)),1)) .^ ratio - 1));
    end

    function par = SS_exact_comp_calculation(x, par)
        alpha_c = par.alpha(x);
        beta_c  = par.beta(x);
        gamma_c = par.gamma(x);

        par.L(par.index.nodes) = -((alpha_c.' * par.Cp) - (beta_c.' * par.Cm));
        par.L(par.index.comps) = -gamma_c.';
    end

    function phi = NLP_obj(x, par)
        index = par.index.pipes(end):par.n.all;
        phi = (1 - param.Econ_weight) * Fval(x) - par.L(index) * x(index);
    end

    function Psi = merit(r, Phi, rho)
        Psi = norm(r, inf) + rho * Phi;
    end

    function Psi_model = model_merit(r_curr, J, s, Phi_curr, par, rho)
        Psi_model = norm(r_curr + J * s, inf) + rho * (Phi_curr + par.L * s);
    end

    function delta = trust_region(param, s, delta, rho_tr)
        expand_factor = 1.25;
        shrink_factor = 0.5;
        min_delta     = 1e-3;
        max_delta     = 5 * param.step_size;

        boundary_ratio = norm(s, inf) / (delta + eps);

        if rho_tr < 0.85
            delta = max(shrink_factor * delta, min_delta);
        elseif rho_tr > 0.85 && boundary_ratio > 0.85
            delta = min(expand_factor * delta, max_delta);
        end
    end
end








