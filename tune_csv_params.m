function [best_param, best_result, all_trials] = tune_csv_params(net, par, param0, x0, varargin)
% TUNE_SLP_PARAMS — step_size tuner for solveSLP
%
% Chooser:
%   1) Compute best_final = min(final feasible objective).
%   2) Keep runs with final feasible objective <= (1+FinalSlack)*best_final.
%   3) Among kept, minimize Score = W_final*rank(final_obj) + W_speed*rank(k_to_near_best) + W_best*rank(best_feas_obj),
%      where k_to_near_best = earliest k with feasible & obj <= (1+SpeedSlack)*best_final.
%   4) Tie-break by fewer total iterations.
%   Fallback A (no kept): among feasible runs, min final_obj → min k_to_near_best → min best_feas_obj.
%   Fallback B (none feasible): min final_violation → min final_obj.
%
% Options:
%   'RefineRounds'      (3)     linear neighbor refinements
%   'PointsPerRound'    (3)     interior points per refinement
%   'MinStep'           (0.01)
%   'MaxStep'           (1000)
%   'Verb'              (true)
%   'EvalMaxIter'       ([] )   force param.max_iter per eval (default 100 if missing)
%   'FeasTolMultiplier' (2)     feasibility gate = mult * 10^(error_tol)
%   'TailLen'           (5)     only for infeasible fallback tie-breaks
%   'StopOnError'       (true)  rethrow solveSLP errors
%   'FinalSlack'        (0.01)  allowed % gap vs best_final (e.g., 1% = 0.01)
%   'SpeedSlack'        (0.02)  threshold to count "near-best" during the run
%   'W_final'           (1.0)   weight for final objective rank
%   'W_speed'           (0.7)   weight for speed-to-near-best rank
%   'W_best'            (0.3)   weight for best feasible objective seen during run

% ---- Parse options ----
p = inputParser;
addParameter(p,'RefineRounds',0);
addParameter(p,'PointsPerRound',10);
addParameter(p,'MinStep',0.01);
addParameter(p,'MaxStep',1000);
addParameter(p,'Verb',true);
addParameter(p,'EvalMaxIter',[]);
addParameter(p,'FeasTolMultiplier',2);
addParameter(p,'TailLen',5);
addParameter(p,'StopOnError',true);
addParameter(p,'FinalSlack',0.01);
addParameter(p,'SpeedSlack',0.02);
addParameter(p,'W_final',1.0);
addParameter(p,'W_speed',1.0);
addParameter(p,'W_best',0.3);
addParameter(p,'MaxBudget',[]); addParameter(p,'R0',[]); addParameter(p,'Eta',[]);
addParameter(p,'NumStepSizes',[]); addParameter(p,'StepMinMax',[]);
addParameter(p,'NomTolGrid',[]); addParameter(p,'ErrTolGrid',[]);
addParameter(p,'KeepFrac',[]);
parse(p,varargin{:}); cfg = p.Results;

% ---- Probe set ----
probe_steps = [0.1 0.5 1 3 5 10 50 100 200 500];
probe_steps = probe_steps(probe_steps >= cfg.MinStep & probe_steps <= cfg.MaxStep);

% Row template
trial_template = struct( ...
    'step_size',[], 'iters',[], ...
    'converged',[], 'feasible',[], 'final_feasible',[], ...
    'best_feas_obj',[], 'k_at_best_feas',[], ...
    'final_objective',[], 'final_violation',[], 'tail_median_violation',[], ...
    'k_to_near_best',[], ...
    'error_id','', 'error_msg','', ...
    'objs',[], 'viols',[], ...
    'x',[], 'history',[], 'par',[]);
trials = repmat(trial_template,0,1);

if cfg.Verb
    fprintf('[tune] Initial probes: '); fprintf('%g ', probe_steps); fprintf('\n');
end

for s = probe_steps
    trials(end+1,1) = tsp_eval_step(s, net, par, param0, x0, cfg); %#ok<AGROW>
end

best_idx = tsp_select_best(trials, cfg);
best = trials(best_idx);
if cfg.Verb, tsp_print_best('After probes', best); end

% ---- Optional neighbor refinement ----
[lo, hi] = tsp_neighbor_bracket(best.step_size, probe_steps, cfg.MinStep, cfg.MaxStep);
for r = 1:cfg.RefineRounds
    if ~(isfinite(lo) && isfinite(hi)) || hi <= lo, break; end
    pts = linspace(lo, hi, max(2+cfg.PointsPerRound,3));
    pts = pts(2:end-1);  % interior only
    pts = pts(~ismembertol(pts, [trials.step_size], 1e-12, 'DataScale', 1));
    if isempty(pts), break; end

    if cfg.Verb
        fprintf('[tune] Refine %d: bracket [%.6g, %.6g] → ', r, lo, hi); fprintf('%g ', pts); fprintf('\n');
    end

    for s = pts
        trials(end+1,1) = tsp_eval_step(s, net, par, param0, x0, cfg); %#ok<AGROW>
    end

    best_idx = tsp_select_best(trials, cfg);
    best = trials(best_idx);

    sampled_now = sort(unique([probe_steps, pts]));
    [lo, hi] = tsp_neighbor_bracket(best.step_size, sampled_now, cfg.MinStep, cfg.MaxStep);

    if cfg.Verb, tsp_print_best(sprintf('Round %d best', r), best); end
end

% ---- Outputs ----
best_param = param0;
best_param.step_size = best.step_size;

% Compute achieved tolerances from best.history
achieved_error_tol = NaN; 
achieved_nom_quality_tol = NaN;
if ~isempty(best.history)
    % Achieved error tol from min constraint_violation across iterations
    if isfield(best.history, 'constraint_violation')
        cv = arrayfun(@(h) double(h.constraint_violation), best.history, 'uniformoutput', true);
        cv = cv(isfinite(cv) & cv > 0);
        if ~isempty(cv)
            min_cv = min(cv);
            achieved_error_tol = ceil(log10(max(min_cv, realmin('double'))));
        end
    end
    % Achieved nom_quality_tol
    if isfield(best.history, 'quality')
        q = arrayfun(@(h) double(h.quality), best.history, 'uniformoutput', true);
        if numel(q) >= 2 && any(isfinite(q))
            dq = abs(diff(q));
            dq = dq(isfinite(dq) & dq > 0);
            if ~isempty(dq)
                min_dq = min(dq);
                achieved_nom_quality_tol = ceil(log10(max(min_dq, realmin('double')))); % e.g., 3e-4 -> -3
            end
        end
    end
end
% Fallbacks if series missing
if ~isfinite(achieved_error_tol),       achieved_error_tol       = param0.error_tol;       end
if ~isfinite(achieved_nom_quality_tol), achieved_nom_quality_tol = param0.nom_quality_tol; end

best_param.error_tol       = achieved_error_tol;
best_param.nom_quality_tol = achieved_nom_quality_tol;

% Return a compact record of the winning run and the trials table
best_result = rmfield(best, intersect({'x','history','par'}, fieldnames(best)));
all_trials  = tsp_to_table(trials);
if cfg.Verb, tsp_print_best('DONE', best); end

% ===================== helpers =====================

function trial = tsp_eval_step(step_size, net, par, param0, x0, cfg)
    param = param0; 
    param.step_size = step_size;
    if ~isfield(param,'error_tol') || isempty(param.error_tol),        param.error_tol = -6; end
    if ~isfield(param,'nom_quality_tol') || isempty(param.nom_quality_tol), param.nom_quality_tol = -6; end
    if ~isfield(param,'max_iter') || isempty(param.max_iter) || param.max_iter < 1
        param.max_iter = 100;
    end
    if ~isempty(cfg.EvalMaxIter), param.max_iter = cfg.EvalMaxIter; end

    % Call solveSLP
    try
        [x, H, par_out] = solveSLP(net, par, param, x0, []);
        err_id=''; err_msg='';
    catch ME
        if cfg.StopOnError
            fprintf(2,'[tune][ERROR] step=%.6g | %s: %s\n', step_size, ME.identifier, ME.message);
            rethrow(ME);
        else
            H=[]; x=[]; par_out=[]; err_id=ME.identifier; err_msg=ME.message;
        end
    end

    if ~isempty(H)
        S = tsp_extract_signals(H, param, cfg);
        trial = struct( ...
            'step_size', step_size, ...
            'iters', S.K, ...
            'converged', S.converged, ...
            'feasible', S.has_any_feas, ...
            'final_feasible', S.final_feasible, ...
            'best_feas_obj', S.best_feas_obj, ...
            'k_at_best_feas', S.k_at_best_feas, ...
            'final_objective', S.final_obj, ...
            'final_violation', S.final_viol, ...
            'tail_median_violation', S.tail_median_viol, ...
            'k_to_near_best', NaN, ... % filled later once best_final is known
            'error_id', err_id, 'error_msg', err_msg, ...
            'objs', S.objs, 'viols', S.viols, ...
            'x', x, 'history', H, 'par', par_out);
    else
        if cfg.Verb, fprintf(2,'[tune][WARN] step=%.6g returned empty history.\n', step_size); end
        trial = struct( ...
            'step_size', step_size, ...
            'iters', 0, ...
            'converged', false, ...
            'feasible', false, ...
            'final_feasible', false, ...
            'best_feas_obj', Inf, ...
            'k_at_best_feas', Inf, ...
            'final_objective', Inf, ...
            'final_violation', Inf, ...
            'tail_median_violation', Inf, ...
            'k_to_near_best', NaN, ...
            'error_id','', 'error_msg','', ...
            'objs', [], 'viols', [], ...
            'x', [], 'history', [], 'par', []);
    end
end

function S = tsp_extract_signals(H, param, cfg)
    K = numel(H);
    viols = nan(K,1); objs = nan(K,1); steps = nan(K,1);
    for k = 1:K
        if isfield(H(k),'constraint_violation') && ~isempty(H(k).constraint_violation)
            v = H(k).constraint_violation; if isfinite(v), viols(k) = v; end
        end
        if isfield(H(k),'total_objective') && ~isempty(H(k).total_objective)
            o = H(k).total_objective; if isfinite(o), objs(k) = o; end
        end
        if isfield(H(k),'step') && ~isempty(H(k).step)
            s = H(k).step; if isfinite(s), steps(k) = s; end
        end
    end

    feas_tol = cfg.FeasTolMultiplier * 10^(param.error_tol);

    feas_set = find(isfinite(viols) & viols < feas_tol & isfinite(objs));
    has_any_feas = ~isempty(feas_set);
    if has_any_feas
        [best_feas_obj, rel] = min(objs(feas_set));
        k_at_best_feas = feas_set(rel);
    else
        best_feas_obj = Inf; k_at_best_feas = Inf;
    end

    last_obj = find(isfinite(objs),1,'last'); final_obj  = iff(isempty(last_obj), Inf, objs(last_obj));
    last_vio = find(isfinite(viols),1,'last'); final_viol = iff(isempty(last_vio), Inf, viols(last_vio));
    last_stp = find(isfinite(steps),1,'last'); final_step = iff(isempty(last_stp), Inf, steps(last_stp));

    final_feasible  = isfinite(final_viol) && final_viol < feas_tol;
    small_step_gate = final_step < max(1e-6, 0.05*param.step_size);
    early_exit      = (K < param.max_iter);
    converged = final_feasible && (early_exit || small_step_gate);

    tailN = max(1, min(cfg.TailLen, K));
    tail = viols(max(1,K-tailN+1):K); tail = tail(isfinite(tail));
    tail_median = iff(isempty(tail), Inf, median(tail));

    S = struct('K',K,'objs',objs,'viols',viols, ...
               'has_any_feas',has_any_feas,'final_feasible',final_feasible,'converged',converged, ...
               'best_feas_obj',best_feas_obj,'k_at_best_feas',k_at_best_feas, ...
               'final_obj',final_obj,'final_viol',final_viol,'tail_median_viol',tail_median);
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end

function [lo, hi] = tsp_neighbor_bracket(s_best, sampled, smin, smax)
    sampled = sort(unique(sampled));
    idx = find(abs(sampled - s_best) <= 1e-12*max(1,abs(s_best)), 1, 'first');
    if isempty(idx), lo=-inf; hi=inf; return; end
    if idx==1
        lo = max(smin, sampled(1));           hi = min(smax, sampled(min(2,end)));
    elseif idx==numel(sampled)
        lo = max(smin, sampled(max(end-1,1))); hi = min(smax, sampled(end));
    else
        lo = max(smin, sampled(idx-1));       hi = min(smax, sampled(idx+1));
    end
    if hi < lo, t=lo; lo=hi; hi=t; end
end

function k = tsp_select_best(T, cfg)
    % Compute best_final among final-feasible runs
    final_feas = arrayfun(@(t)t.final_feasible, T);
    final_objs = arrayfun(@(t)t.final_objective, T);
    if any(final_feas)
        best_final = min(final_objs(final_feas));
        % k_to_near_best: first k with feas & obj <= (1+SpeedSlack)*best_final
        near_thresh = (1 + cfg.SpeedSlack) * best_final;
        feas_tol = [];
        for i = 1:numel(T)
            if isempty(T(i).objs), T(i).k_to_near_best = Inf; continue; end
            objs = T(i).objs; viols = T(i).viols;
            feas_gate = cfg.FeasTolMultiplier * 10^(-6); % fallback gate; param.error_tol not stored per trial
            feas_mask = isfinite(viols) & (viols < feas_gate);
            near_mask = isfinite(objs) & (objs <= near_thresh);
            idx = find(feas_mask & near_mask, 1, 'first');
            if isempty(idx), T(i).k_to_near_best = Inf; else, T(i).k_to_near_best = idx; end
        end
        % Keep runs within FinalSlack of best_final
        keep = final_feas & (final_objs <= (1 + cfg.FinalSlack) * best_final);
        if any(keep)
            k = rank_composite(T(keep), cfg);
            keep_idx = find(keep);
            k = keep_idx(k);
            return
        else
            % Fallback A: among feasible runs, prefer lowest final objective -> fastest to near-best -> best_feas_obj
            feas_any = arrayfun(@(t)t.feasible, T);
            if any(feas_any)
                idx = find(feas_any);
                k = rank_feasible_fallback(T(idx), cfg);
                k = idx(k);
                return
            end
        end
    end

    % Fallback B: none feasible — min final violation -> min final objective
    vfin = arrayfun(@(t)t.final_violation, T);
    [~, r1] = sort(vfin, 'ascend');
    [~, r2] = sort(final_objs, 'ascend');
    score = rankpos(r1) + 0.25*rankpos(r2);
    [~, k] = min(score);
end

function k_rel = rank_composite(S, cfg)
    fobj = arrayfun(@(t)t.final_objective, S);
    k2nb = arrayfun(@(t)t.k_to_near_best, S);
    bobj = arrayfun(@(t)t.best_feas_obj, S);
    itrs = arrayfun(@(t)t.iters, S);
    [~, rF]  = sort(fobj, 'ascend');
    [~, rS]  = sort(k2nb, 'ascend');
    [~, rB]  = sort(bobj, 'ascend');
    score = cfg.W_final*rankpos(rF) + cfg.W_speed*rankpos(rS) + cfg.W_best*rankpos(rB);
    [~, rI]  = sort(itrs, 'ascend'); score = score + 0.1*rankpos(rI);
    [~, k_rel] = min(score);
end

function k_rel = rank_feasible_fallback(S, cfg)
    fobj = arrayfun(@(t)t.final_objective, S);
    bobj = arrayfun(@(t)t.best_feas_obj, S);
    itrs = arrayfun(@(t)t.iters, S);
    k2nb = arrayfun(@(t)t.k_to_near_best, S);
    [~, rF] = sort(fobj, 'ascend');
    [~, rS] = sort(k2nb, 'ascend');
    [~, rB] = sort(bobj, 'ascend');
    score = rankpos(rF) + 0.5*rankpos(rS) + 0.25*rankpos(rB);
    [~, rI] = sort(itrs, 'ascend'); score = score + 0.1*rankpos(rI);
    [~, k_rel] = min(score);
end

function rp = rankpos(order_idx)
    N = numel(order_idx); rp = zeros(N,1); rp(order_idx) = 1:N;
end

function T = tsp_to_table(tr)
    if isempty(tr), T = table(); return; end
    n = numel(tr);
    step=zeros(n,1); it=zeros(n,1); conv=false(n,1); feas=false(n,1); ffeas=false(n,1);
    bestobj=zeros(n,1); kbest=zeros(n,1); objfin=zeros(n,1);
    vfin=zeros(n,1); vmed=zeros(n,1); k2nb=zeros(n,1);
    for i=1:n
        step(i)=tr(i).step_size; it(i)=tr(i).iters;
        conv(i)=logical(tr(i).converged); feas(i)=logical(tr(i).feasible); ffeas(i)=logical(tr(i).final_feasible);
        bestobj(i)=tr(i).best_feas_obj; kbest(i)=tr(i).k_at_best_feas; objfin(i)=tr(i).final_objective;
        vfin(i)=tr(i).final_violation; vmed(i)=tr(i).tail_median_violation;
        k2nb(i)=tr(i).k_to_near_best;
    end
    T = table(step,it,conv,feas,ffeas,bestobj,kbest,objfin,vfin,vmed,k2nb, ...
        'VariableNames',{'step_size','iters','converged','feasible','final_feasible','best_feas_obj','k_at_best_feas','final_objective','final_violation','tail_median_violation','k_to_near_best'});
end

function tsp_print_best(tag, b)
    fprintf('[tune] %s: step=%.6g | CONV=%d | FINAL_FEAS=%d | final_obj=%.6g | k_to_near_best=%s | best_feas_obj=%.6g | iters=%d\n', ...
        tag, b.step_size, b.converged, b.final_feasible, b.final_objective, tsp_kstr(b.k_to_near_best), b.best_feas_obj, b.iters);
end

function s = tsp_kstr(k), if isfinite(k), s=sprintf('%d',k); else, s='—'; end, end

end
