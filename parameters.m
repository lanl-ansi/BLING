function [par, x, history] = parameters(net, bc, param)
% SS_PARAMETERS  Build problem structures for the SLP/Newton pipeline.
% Returns:
%   par     : all matrices/indices
%   x       : initial guess
%   history : struct array with fields:
%             F_comp_val, F_surp_val, delta, step, constraint_violation,
%             lambda_eq, lambda_ineq, duality_gap, quality

par = struct();

% ----------------------------- Indices
par = indexing(par);

% ----------------------------- Incidence & resistance
[par, Cmax, Cmin] = incidence(par);

% ----------------------------- BCs, scaling & bounds
[par, slackPressure] = bc_scale_bounds(par, Cmax, Cmin);

% ----------------------------- Linear objective
par = objective(par);

% ----------------------------- Linear Jacobian parts
[par, Abase, Anewton]  = linear_jacobian(par, slackPressure);

% ----------------------------- Nonlinear Jacobian parts
par = nonlinear_jacobian(par, Abase, Anewton);

% ----------------------------- Initial guess
x = initial_guess(par);

% ----------------------------- Preallocate history
template = struct( ...
    'F_comp_val', NaN, ...
    'F_surp_val', NaN, ...
    'delta', NaN, ...
    'step', NaN, ...
    'constraint_violation', NaN, ...
    'lambda_gnode', [], ...
    'lambda_all', [], ...
    'duality_gap', NaN, ...
    'quality', NaN, ...
    'total_objective', NaN, ...
    'comp_time',NaN);

history = repmat(template, param.max_iter, 1);

% ========================= Nested functions =========================

% ---- Indices
    function par = indexing(par)
        % Sizes
        par.n.nodes  = size(net.nodes,  1);
        par.n.comps  = size(net.comps,  1);
        par.n.pipes  = size(net.pipes,  1);
        par.n.gnodes = size(net.gnodes, 1);

        % Decision count: [nodes | comps | pipes | gnode_s | gnode_d]
        par.n.all = par.n.nodes + par.n.comps + par.n.pipes + 2*par.n.gnodes;

        % Variable indices (column vectors)
        par.index.nodes   = (1:par.n.nodes).';
        par.index.comps   = par.n.nodes + (1:par.n.comps).';
        par.index.pipes   = par.n.nodes + par.n.comps + (1:par.n.pipes).';
        par.index.gnode_s = par.n.nodes + par.n.comps + par.n.pipes + (1:par.n.gnodes).';
        par.index.gnode_d = par.n.nodes + par.n.comps + par.n.pipes + par.n.gnodes + (1:par.n.gnodes).';

        % Slack node indices
        par.index.slack = find(net.nodes(:,end));
        par.n.slack     = numel(par.index.slack);

        % Newton variable count (nodes + comps + pipes)
        par.n.newton = par.n.nodes + par.n.comps + par.n.pipes;
    end

% ---- Incidence & resistance
    function [par, Cmax, Cmin] = incidence(par)
        % Pipes
        p_id  = net.pipes(:,1);
        p_fr  = net.pipes(:,2);
        p_to  = net.pipes(:,3);

        par.P = sparse([p_fr;  p_to], [p_id; p_id], ...
            [-ones(par.n.pipes,1); ones(par.n.pipes,1)], ...
            par.n.nodes, par.n.pipes).';

        % Compressors
        c_id  = net.comps(:,1);
        c_fr  = net.comps(:,2);
        c_to  = net.comps(:,3);
        c_min = net.comps(:,4); 
        c_max = net.comps(:,5);  

        par.C = sparse([c_fr;  c_to], [c_id;  c_id], ...
            [-ones(par.n.comps,1);  ones(par.n.comps,1)], ...
            par.n.nodes, par.n.comps).';

        Cmax  = sparse([c_fr;  c_to], [c_id;  c_id], ...
            [-c_max;  ones(par.n.comps,1)], ...
            par.n.nodes, par.n.comps).';

        Cmin  = sparse([c_fr;  c_to], [c_id;  c_id], ...
            [-c_min;  ones(par.n.comps,1)], ...
            par.n.nodes, par.n.comps).';

        % Positive/negative selector for compressors
        par.Cp = (par.C > 0);
        par.Cm = (par.C < 0);

        % G-node selector: G : [nNodes x nGnodes], 1 if gnode j mapped to node i
        par.G = sparse(net.gnodes(:,2), net.gnodes(:,1), 1, par.n.nodes, par.n.gnodes);

        % Pipe hydraulic resistance Lambda
        diam   = net.pipes(:,4);
        area   = (pi/4) * diam.^2;
        length = net.pipes(:,5);
        fric   = net.pipes(:,6);

        par.Lambda = fric .* length ./ (2 * diam .* area.^2);
        par.Lambda = par.Lambda(:);
    end

% ---- Boundary conditions, scaling, and bounds
    function [par, slackPressure]  = bc_scale_bounds(par, Cmax, Cmin)
        % Slack selector S : [nSlack x nNodes]
        par.S = sparse(1:par.n.slack, par.index.slack, 1, par.n.slack, par.n.nodes);
        slackPressure = bc.pslack(2,:).';

        % Withdrawals/injections at nodes
        g = par.G * bc.gbar(2,:).';
        par.q = g + bc.qbar(2,:).';

        % Bounds from data
        pmin = net.nodes(:,4);
        pmax = net.nodes(:,5);
        fcmin = net.comps(:,6);
        fcmax = net.comps(:,7);
        smax = bc.smax(2,:).';
        dmax = bc.dmax(2,:).';
        smin = zeros(size(smax));
        dmin = zeros(size(dmax));

        % ---- Nondimensional scaling
        universalR = 8314.472; molecwghtAir = 28.9626; gasG = param.Gas_specific_gravity;

        R = universalR / (molecwghtAir * gasG);
        Z = param.compressibility_factor;
        T = param.Temperature;

        par.sigma = sqrt(max(Z,eps) * max(R,eps) * max(T,eps));

        % Use a robust xi tied to Lambda
        Lmean = mean(par.Lambda);
        Lmax = (max(par.Lambda)+Lmean)/4;
        par.xi = 1 / sqrt(max(Lmax, eps));

        % Pressure scaling
        ps_ref = min(slackPressure);
        par.varp = ps_ref;
        par.varq = par.varp * par.xi / par.sigma;

        % Scale Lambda to nondimensional
        par.Lambda = par.Lambda * par.xi^2;

        % Scale everything else
        slackPressure = slackPressure / par.varp;
        par.q = par.q / par.varq;
        g     = g     / par.varq;

        pmax = pmax / par.varp;    pmin = pmin / par.varp;
        smax = smax / par.varq;    dmax = dmax / par.varq;
        fcmin = fcmin / par.varq;  fcmax = fcmax / par.varq;

        % ---- Decision variable bounds
        par.xmin = -Inf(par.n.all, 1);
        par.xmax =  Inf(par.n.all, 1);
        par.pmin = pmin;

        Pp = (par.P > 0);   % inflow node for each pipe
        Pm = (par.P < 0);   % outflow node for each pipe

        % nodes
        par.xmin(par.index.nodes) = pmin;
        par.xmax(par.index.nodes) = pmax;

        pmin2 = pmin.^2;
        pmax2 = pmax.^2;
        num = 0.5 * (Pp * pmax2 - Pm * pmin2); 
        num = max(num, 0);         
        den = max(par.Lambda, eps);  
        flowcap = sqrt(num ./ den);  
        par.xmin(par.index.pipes) = - flowcap;
        par.xmax(par.index.pipes) =   flowcap;

                % compressors
        par.xmin(par.index.comps) = fcmin;
        par.xmax(par.index.comps) = min(max(flowcap), fcmax);

        % g-nodes
        par.xmin([par.index.gnode_s; par.index.gnode_d]) = [smin; dmin];
        par.xmax([par.index.gnode_s; par.index.gnode_d]) = [smax; dmax];

        % ---- Linear inequality constraints for compressors
        par.B = spalloc(2*par.n.comps, par.n.all, 4*par.n.comps);
        par.B(:, par.index.nodes) = [-Cmin; Cmax];
        par.b = zeros(2*par.n.comps, 1);
    end

% ---- Global linear objective (overwritten by SLP if exact_comp_energy==1)
    function par = objective(par)
        r = param.specific_heat_capacity_ratio;
        ratio = (r - 1) / r;

        xPmax = par.xmax(par.index.nodes);
        xPmin = par.xmin(par.index.nodes);
        denom = max(xPmin, eps);

        alpha = ((1 - r) / r) ...
            * mean(par.xmax(par.index.comps)) / max(mean(xPmax), eps) ...
            * min(xPmax ./ denom).^ratio;

        gamma = ( mean(xPmax ./ denom).^ratio - 1 );

        Lp  = (1 - param.Econ_weight) * alpha * sum(par.C, 1);
        Lfc = (1 - param.Econ_weight) * gamma * ones(1, par.n.comps);
        Ls  =  param.Econ_weight * bc.cs(2,:);
        Ld  =  param.Econ_weight * bc.cd(2,:);

        par.L = [-Lp, Lfc, zeros(1,par.n.pipes), Ls, -Ld];
    end

% ---- Linear matrix elements of Jacobians
    function [par, Abase, Anewton] = linear_jacobian(par, slackPressure)
        nRows = par.n.slack + par.n.pipes + par.n.nodes;
        nCols = par.n.nodes + par.n.comps + par.n.pipes + 2*par.n.gnodes;
        nVals = par.n.slack + 4*par.n.pipes + 2*par.n.gnodes;

        if param.newton_solve == 1
            nCols_newton = par.n.nodes + par.n.comps + par.n.pipes;
            nRows_newton = par.n.nodes + par.n.comps + par.n.pipes + par.n.slack;
            nVals_newton = par.n.slack + 4*par.n.pipes + par.n.slack;
        end

        % Initialize linear matrix terms
        Abase = spalloc(nRows, nCols, nVals);
        par.a = zeros(nRows, 1);

        if param.newton_solve == 1
            Anewton = spalloc(nRows_newton, nCols_newton, nVals_newton);
            par.index.newton = [par.index.nodes; par.index.comps; par.index.pipes];
        else
            Anewton = [];
        end

        % ---- Fill Jacobian matrices
        % Slack pressure rows
        row = 1:par.n.slack;
        Abase(row, par.index.nodes) = par.S;
        if param.newton_solve == 1
            Anewton(row, par.index.nodes) = par.S;
        end
        par.a(row) = slackPressure;

        % Nodal balance rows
        row = par.n.slack + par.n.pipes + (1:par.n.nodes);
        Abase(row, [par.index.comps; par.index.pipes; par.index.gnode_s; par.index.gnode_d]) = ...
            [par.C.', par.P.', par.G, -par.G];

        if param.newton_solve == 1
            Anewton(row, [par.index.comps; par.index.pipes]) = [par.C.', par.P.'];
        end

        par.a(row) = par.q;
    end

% ---- Nonlinear matrix elements of Jacobians
    function par = nonlinear_jacobian(par, Abase, Anewton)
        nRows    = par.n.slack + par.n.pipes + par.n.nodes;
        pipeRows = par.n.slack + (1:par.n.pipes);
        pipeCols = par.index.pipes;

        g = 9.80665;
        len  = net.pipes(:,5);
        sina = sind(net.pipes(:,end));
        coef = g / (2 * par.sigma^2);
        incline = coef * [len .* sina; len .* sina];

        % EOS coefficients
        if get_field(param,'ideal_eos',1)==0
            p_atm = 101350;  Gs = 0.650784;
            a1 = 344400; a2 = 1.785; a3 = 3.825;
            T  = param.Temperature;
            powT = (max(1.8*T, eps))^a3;
            A    = a1 * 10^(a2 * Gs) / powT;
            inv689 = 1/6894.75729;
            tmp = (p_atm*inv689) * A;
            b1 = 1 + tmp;
            b2 = (par.varp*inv689) * A;
        else
            b1 = 1; b2 = 0;
        end
        b1 = b1 / 2; b2 = b2 / 3;

        % For sparse add-ins
        PipeRow = par.n.slack + [net.pipes(:,1); net.pipes(:,1)];
        PipeCol = [net.pipes(:,2); net.pipes(:,3)]; 

        % Node pressure nonlinear terms
        in  = par.index.nodes(net.pipes(:,2));
        out = par.index.nodes(net.pipes(:,3));
        p1 = @(x,a) a * [-x(in);  x(out)];
        p2 = @(x,b) b * [-x(in).^2; x(out).^2];
        p3 = @(x,c) c * [x(in).^3; x(out).^3];
        p_drop  = @(x,a,b) b1 .* p1(x,a) + b2 .* p2(x,b);
        p_angle = @(x,a,b,c) incline .* ( b1.^2 .* abs(p1(x,a)) + 2 * b1 .* b2 .* abs(p2(x,b)) + b2.^2 .* p3(x,c) );

        % Equality constraints and Jacobian of NLP
        par.A = @(x,a,b,c) Abase + ...
            sparse(PipeRow, PipeCol, p_drop(x,a,b) - p_angle(x,a,b,c), nRows, par.n.all) + ...
            sparse(pipeRows, pipeCols, par.Lambda .* (abs(a) * abs(x(pipeCols))), nRows, par.n.all);

        % Equations and Jacobian for Newton solve
        if get_field(param,'newton_solve',0) == 1
            nRows_newton = par.n.nodes + par.n.comps + par.n.pipes + par.n.slack;
            par.A_newton = @(x,a,b,c) Anewton + ...
                sparse(PipeRow, PipeCol, p_drop(x,a,b) - p_angle(x,a,b,c), nRows_newton, par.n.newton) + ...
                sparse(pipeRows, pipeCols, par.Lambda .* (abs(a) * abs(x(pipeCols))), nRows_newton, par.n.newton);
        end
    end

% ---- Initial guess
    function x = initial_guess(par)
        x = 0.05 + 0.03 * rand(par.n.all, 1);
        x(par.index.nodes)   = 0.5 * (par.xmax(par.index.nodes) + par.xmin(par.index.nodes));
        x(par.index.gnode_d) = 0.5 * par.xmax(par.index.gnode_d);
    end

% ---- Utility
    function v = get_field(S, name, default)
        if isfield(S,name) && ~isempty(S.(name)), v = S.(name);
        else, v = default; end
    end
end