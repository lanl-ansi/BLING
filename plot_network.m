function h = plot_network(net, par, param, lambda, varargin)
% PLOT_GAS_NETWORK  Gas network plot (flows on edges, pressures on nodes, and prices at nodes).
%   H = PLOT_GAS_NETWORK(NET, PAR, PARAM, LAMBDA)
%   H = PLOT_GAS_NETWORK(NET, PAR, PARAM, LAMBDA, 'Name',Value,...)
%   H = PLOT_GAS_NETWORK(NET, PAR, PARAM, LAMBDA, X, 'Name',Value,...)
%
% Behavior:
%   - If X is omitted   -> preview mode (topology-only, very fast), flags ignored
%   - If X is provided  -> up to FIVE separate figures:
%       * param.solution_plot == true    -> "solution"      (flows + pressures)
%       * param.pressure_plot == true    -> "pressure"      (nodes by P; black edges)
%       * param.flow_plot/.plot_flow     -> "flow"          (edges by |m|; G-nodes by |q+Gxd-Gxs|)
%       * param.gnode_flow_plot == true  -> "gnode_flow"    (edges black; subset of nodes colored)
%       * param.price_plot == true       -> "price"         (edges black; nodes by price from lambda)
%
% Name-Value Options
%   'flow_width'        : [wmin wmax] (default [0.5 2.5])
%   'phi_target'        : node packing fraction (default 0.06)
%   'min_sep'           : min node separation (default 0.06 of span)
%   'iter_triOL'        : compressor triangle de-overlap iters (default 8)
%   'tri_scale'         : global triangle side multiplier (default 1.50)
%   'tri_small_q'       : quantile for "small nodes" mean (default 0.25)
%   'AxesPosition'      : axes position (default [0.12 0.20 0.70 0.72])
%   'gnode_min_factor'  : (<1) G-node min radius vs triangle min circumR (default 0.6)
%   'ModelName'         : optional title string
%   'num_width_bins'    : edge-width bin count for batching (default 5)

%% ---------- accept optional X then name/value ----------
x = [];
if ~isempty(varargin) && (isnumeric(varargin{1}) || islogical(varargin{1}))
    x = varargin{1};
    varargin(1) = [];
end
 
p = inputParser; p.CaseSensitive = false;
addParameter(p,'flow_width',[0.5 2.5]);
addParameter(p,'phi_target',0.06);
addParameter(p,'min_sep',0.06);
addParameter(p,'iter_triOL',8);
addParameter(p,'tri_scale',1.50);
addParameter(p,'tri_small_q',0.25);
addParameter(p,'AxesPosition',[0.12 0.20 0.70 0.72]);
addParameter(p,'gnode_min_factor',0.6);
addParameter(p,'ModelName','');
addParameter(p,'num_width_bins',5);
parse(p, varargin{:});
opt = p.Results;

%% ---------- which figures to make ----------
hasX = ~isempty(x);
flag_solution = isfield(param,'pressure_flow_plot') && param.pressure_flow_plot;
flag_pressure = isfield(param,'pressure_plot')      && param.pressure_plot;
flag_flow     = (isfield(param,'flow_plot') && param.flow_plot) || ...
    (isfield(param,'plot_flow') && param.plot_flow);
flag_gnflow   = isfield(param,'gnode_flow_plot')    && param.gnode_flow_plot;
flag_price    = isfield(param,'LMP_plot')           && param.LMP_plot;

modes = {};
if ~hasX
    modes = {'topology'};
else
    if flag_solution, modes{end+1} = 'solution';   end
    if flag_pressure, modes{end+1} = 'pressure';   end
    if flag_flow,     modes{end+1} = 'flow';       end
    if flag_gnflow,   modes{end+1} = 'gnode_flow'; end
    if flag_price,    modes{end+1} = 'price';      end
    if isempty(modes), modes = {'solution'}; end
end

% ----------- FONT SIZE PRESET -----------
FS.axes    = 28;
FS.title   = 32;
FS.cbTicks = 20;
FS.cbLabel = 28;

% Render one figure per requested mode
h = struct();
for mi = 1:numel(modes)
    h.(modes{mi}) = render_one(net, par, x, lambda, opt, modes{mi}, FS);
end
end % ======= end main =======


% ==============================================================
% =============== internal renderer (one figure) ================
% ==============================================================

function h = render_one(net, par, x, lambda, opt, mode, FS)
% mode in {'topology','solution','pressure','flow','gnode_flow','price'}

%% ---------- large number of nodes threshold ----------
BIG_N  = 1e5;

%% ---------- extract network ----------
Xc = net.nodes(:,2);  Yc = net.nodes(:,3);  N = numel(Xc);
pipes = net.pipes(:,[2,3]);    nP = size(pipes,1);
comps = net.comps(:,[2,3]);    nC = size(comps,1);

Gext = [];
if isfield(net,'gnodes') && ~isempty(net.gnodes), Gext = net.gnodes(:,2); end
slack_ids = par.index.slack(:);

% external id -> compact index
ext_ids = net.nodes(:,1); maxid = max(ext_ids);
map  = sparse(ext_ids,1,1:N,maxid,1);
pf   = full(map(pipes(:,1))); pt = full(map(pipes(:,2)));
cf   = full(map(comps(:,1))); ct = full(map(comps(:,2)));
Gi   = full(map(Gext)); Gi = Gi(Gi>0 & Gi<=N);
slack_idx = full(map(slack_ids));
slack_idx = slack_idx(slack_idx>0);

%% ---------- adapt options for big cases ----------
num_width_bins0 = opt.num_width_bins;
if N > BIG_N
    opt.num_width_bins = min(num_width_bins0, 3);
end

%% ---------- pull data from x (if given) ----------
flowP = []; flowC = []; pressure = [];
gs_vec = []; gd_vec = [];
gnode_net_at_nodes = zeros(N,1);

hasX = ~isempty(x);
if hasX
    if isfield(par.index,'pipes')   && ~isempty(par.index.pipes),   flowP    = x(par.index.pipes);   end
    if isfield(par.index,'comps')   && ~isempty(par.index.comps),   flowC    = x(par.index.comps);   end
    if isfield(par.index,'nodes')   && ~isempty(par.index.nodes),   pressure = x(par.index.nodes);   end

    if isfield(par,'G') && ~isempty(par.G)
        if isfield(par.index,'gnode_s') && ~isempty(par.index.gnode_s)
            gs_vec = par.G * x(par.index.gnode_s);
        end
        if isfield(par.index,'gnode_d') && ~isempty(par.index.gnode_d)
            gd_vec = par.G * x(par.index.gnode_d);
        end
    end

    % ---------- net = q + G*xd - G*xs ----------
    net_ext = zeros(maxid,1);
    if isfield(par,'q') && ~isempty(par.q)
        qv = par.q(:);
        if numel(qv)==N
            tmp_ext = zeros(maxid,1); tmp_ext(ext_ids) = qv;
            net_ext = net_ext + tmp_ext;
        elseif numel(qv)==maxid
            net_ext = net_ext + qv;
        end
    end
    if isfield(par,'G') && ~isempty(par.G)
        if isfield(par.index,'gnode_d') && ~isempty(par.index.gnode_d)
            gd_ext = par.G * x(par.index.gnode_d);
            gd_ext = pad_or_pass_ext(gd_ext, maxid, ext_ids);
            net_ext = net_ext + gd_ext;
        end
        if isfield(par.index,'gnode_s') && ~isempty(par.index.gnode_s)
            gs_ext = par.G * x(par.index.gnode_s);
            gs_ext = pad_or_pass_ext(gs_ext, maxid, ext_ids);
            net_ext = net_ext - gs_ext;
        end
    end
    net_comp = zeros(N,1);
    nz  = find(net_ext ~= 0);
    ci  = full(map(nz));
    keep = (ci > 0 & ci <= N);
    ci  = ci(keep);  nv = net_ext(nz(keep));
    if ~isempty(ci)
        [uc,~,ic] = unique(ci);
        net_comp(uc) = accumarray(ic, nv, [], @sum, 0);
    end
    gnode_net_at_nodes = abs(net_comp);
end
isPreview = strcmp(mode,'topology');

%% ---------- Node prices directly from lambda ----------
price_ext = [];
if ~isempty(lambda)
    lam = lambda(:);
    if numel(lam) == maxid
        price_ext = lam;
    elseif numel(lam) == N
        price_ext = zeros(maxid,1);
        price_ext(ext_ids) = lam;
    else
        price_ext = zeros(maxid,1);
        n = min(numel(lam), maxid);
        price_ext(1:n) = lam(1:n);
    end
else
    price_ext = zeros(maxid,1);
end
price_compact = price_ext(ext_ids);

%% ---------- palettes ----------
if hasX
    if N > BIG_N
        cmapF = cmapFlowTurquoiseBlue(128);
        cmapP = cmapPressureYellowRed(128);
    else
        cmapF = cmapFlowTurquoiseBlue(256);
        cmapP = cmapPressureYellowRed(256);
    end
else
    if N > BIG_N
        cmapF = turbo(128);
        cmapP = parula(128);
    else
        cmapF = turbo(256);
        cmapP = parula(256);
    end
end
colEdgePipe_default = [0.2 0.2 0.2];
colEdgeComp_default = [0.67 0.13 1.0];
colEdgeNode_default = [0.5 0.5 0.5];
colEdgeG            = [0 0.99 1.0];
nodeFill            = [1 1 1];

if strcmp(mode,'pressure') || strcmp(mode,'gnode_flow') || strcmp(mode,'price')
    colEdgePipe = [0.8 0.8 0.8];
    colEdgeComp = [0.8 0.8 0.8];
else
    colEdgePipe = colEdgePipe_default;
    colEdgeComp = colEdgeComp_default;
end

%% ---------- node radii (adaptive) ----------
use_fast_radii = (N > 1e5) && (~strcmp(mode,'flow'));
if use_fast_radii
    rx=max(Xc)-min(Xc); ry=max(Yc)-min(Yc); boxA=max(rx*ry,eps);
    nodeR = 0.5*sqrt((opt.phi_target*boxA)/(N*pi))*ones(N,1);
else
    nodeR = node_radii_knn_fast(Xc, Yc, 1);
end

%% ---------- compressor triangles (adaptive) ----------
simpleTris = (N > BIG_N) || (nC > 2e4);
[triCenters, triSide, triMinSide] = compressor_triangles_adaptive_fast( ...
    Xc, Yc, cf, ct, nodeR, Gi, opt.iter_triOL, opt.tri_scale, opt.tri_small_q, simpleTris);

if ~isempty(Gi) && ~isempty(triMinSide) && all(isfinite(triMinSide))
    Rtri_min = triMinSide / sqrt(3);
    rG_min   = opt.gnode_min_factor * Rtri_min;
    nodeR(Gi) = max(nodeR(Gi), rG_min);
end

%% ---------- axes ----------
figure('Position',[1,94,1128,772],'Renderer','opengl','GraphicsSmoothing','off');
clf;
ax = axes('Units','normalized','Position',opt.AxesPosition);
if ~isempty(opt.ModelName)
    title(ax, opt.ModelName, 'Interpreter','none', 'FontSize', FS.title);
end
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
h = struct();

%% ---------- flow mapping ----------
doFlow_raw = (~isempty(flowP) && numel(flowP)==nP) || ...
    (~isempty(flowC) && numel(flowC)==nC) || ...
    (hasX && any(gnode_net_at_nodes~=0)) || ...
    (~isempty(gs_vec) || ~isempty(gd_vec));

if isPreview
    doFlow = false;
elseif strcmp(mode,'pressure') || strcmp(mode,'gnode_flow') || strcmp(mode,'price')
    doFlow = false;
else
    doFlow = doFlow_raw;
end

%% ---------- unified flow clim ----------
if doFlow
    vals = [];
    if ~isempty(flowP), vals = [vals; abs(flowP(:))]; end
    if ~isempty(flowC), vals = [vals; abs(flowC(:))]; end
    if ~isempty(Gi) && any(gnode_net_at_nodes(Gi) ~= 0)
        vals = [vals; gnode_net_at_nodes(Gi)];
    end
    pr = prctile_safe(vals,[5 95]);
    climF = [0 max(pr(2), max(vals)*1e-6)];
    fw = @(f) opt.flow_width(1) + (opt.flow_width(2)-opt.flow_width(1)) * (min(abs(f),climF(2))/climF(2));
else
    fw = @(f) opt.flow_width(1);
    climF = [];
end

%% ---------- PIPES ----------
if nP>0
    X1p = Xc(pf);  Y1p = Yc(pf);
    X2p = Xc(pt);  Y2p = Yc(pt);

    if isPreview
        Cuni = repmat(colEdgePipe, nP, 1);
        h.pipes = draw_segments_patch(ax, X1p,Y1p,X2p,Y2p, Cuni, opt.flow_width(1));
    elseif doFlow && ~isempty(flowP) && numel(flowP)==nP
        wP = fw(flowP(:));
        cP = map_to_cmap(abs(flowP(:)), climF, cmapF);
        if nP > BIG_N
            h.pipes = draw_segments_binned(ax, X1p, Y1p, X2p, Y2p, abs(flowP(:)), wP, 12, opt.num_width_bins, cmapF);
        else
            [h.pipes,~] = draw_grouped_segments(ax, X1p, Y1p, X2p, Y2p, cP, wP, opt.num_width_bins);
        end
    else
        [XS,YS] = seg_nan([X1p X2p], [Y1p Y2p]);
        h.pipes = plot(ax,XS,YS,'-','Color',colEdgePipe,'LineWidth',opt.flow_width(1),'HitTest','off');
    end
else
    h.pipes = gobjects(0,1);
end

%% ---------- COMPRESSOR EDGES ----------
if nC>0
    X1c = Xc(cf);  Y1c = Yc(cf);
    X2c = Xc(ct);  Y2c = Yc(ct);

    if isPreview
        Cuni = repmat(colEdgePipe, nC, 1);
        h.comp_edges = draw_segments_patch(ax, X1c, Y1c, X2c, Y2c, Cuni, opt.flow_width(1));
    elseif doFlow && ~isempty(flowC) && numel(flowC)==nC
        wC = fw(flowC(:));
        cC = map_to_cmap(abs(flowC(:)), climF, cmapF);
        if nC > BIG_N
            h.comp_edges = draw_segments_binned(ax, X1c, Y1c, X2c, Y2c, abs(flowC(:)), wC, 12, opt.num_width_bins, cmapF);
        else
            [h.comp_edges,~] = draw_grouped_segments(ax, X1c, Y1c, X2c, Y2c, cC, wC, opt.num_width_bins);
        end
    else
        [XS,YS] = seg_nan([X1c X2c], [Y1c Y2c]);
        h.comp_edges = plot(ax,XS,YS,'-','Color',colEdgePipe,'LineWidth',opt.flow_width(1),'HitTest','off');
    end
else
    h.comp_edges = gobjects(0,1);
end

%% ---------- node rendering by mode ----------
useScatterNodes = false; % (N > BIG_N);

if strcmp(mode,'flow') && hasX
    % FLOW MODE (mass flow plot):
    % - Pipes/comp edges already colored by |m|
    % - Plot ONLY G-nodes, colored by | q + G*xd - G*xs | using SAME cmap/limits as pipes

    % Keep legacy fields, but do not draw "all nodes"
    h.nodes_fill   = gobjects(1);
    h.nodes_edge   = gobjects(1);

    if ~isempty(Gi)
        gvals = gnode_net_at_nodes(Gi);           % | q + G*xd - G*xs | at G-nodes
        gRGB  = map_to_cmap(gvals, climF, cmapF);

        if useScatterNodes
            Sg = dataR_to_points2(ax, nodeR(Gi), Xc(Gi), Yc(Gi));

            % white fill for G-nodes only
            h.nodes_fill_G = scatter(ax, Xc(Gi), Yc(Gi), Sg, [1 1 1], ...
                'filled','Marker','o','MarkerEdgeColor','none','HitTest','off');

            % colored border for G-nodes only
            h.nodes_edge_G = scatter(ax, Xc(Gi), Yc(Gi), Sg*1.10, gRGB, ...
                'o','MarkerFaceColor','none','MarkerEdgeColor','flat', ...
                'LineWidth',0.75,'HitTest','off');
        else
            [h.nodes_fill_G, h.nodes_edge_G] = draw_node_circles_batched( ...
                ax, Xc(Gi), Yc(Gi), nodeR(Gi), gRGB, [1 1 1]);
        end
    else
        h.nodes_fill_G = gobjects(1);
        h.nodes_edge_G = gobjects(1);
    end

elseif strcmp(mode,'pressure') && hasX && ~isempty(pressure) && numel(pressure)==N
    prP = prctile_safe(pressure(:),[5 95]); climP = [prP(1) prP(2)];
    if ~(isfinite(climP(1)) && isfinite(climP(2))) || climP(2)<=climP(1)
        climP = [min(pressure) max(pressure)+eps];
    end
    nodeEdgeColors = map_to_cmap(pressure(:),climP,cmapP);
    if useScatterNodes
        S = dataR_to_points2(ax, nodeR, Xc, Yc);
        h.nodes_fill = scatter(ax, Xc, Yc, S, [1 1 1], 'filled', 'Marker','o', 'MarkerEdgeColor','none','HitTest','off');
        h.nodes_edge = scatter(ax, Xc, Yc, S*1.10, nodeEdgeColors, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor','flat','LineWidth',0.75,'HitTest','off');
    else
        [h.nodes_fill, h.nodes_edge] = draw_node_circles_batched(ax, Xc, Yc, nodeR, nodeEdgeColors, [1 1 1]);
    end
    setappdata(ax,'climP',climP);

elseif strcmp(mode,'gnode_flow') && hasX
    sVals = []; dVals = []; sIdx = []; dIdx = [];
    Rmin = 0.175 * max(nodeR);
    if ~isempty(gs_vec), sVals = abs(gs_vec(:)); sIdx = find(sVals>0); Rs = max(nodeR(sIdx), Rmin); end
    if ~isempty(gd_vec), dVals = abs(gd_vec(:)); dIdx = find(dVals>0); Rd = max(nodeR(dIdx), Rmin); end

    if ~isempty(sIdx)
        sClim = local_clim(sVals(sIdx));
        sRGB  = map_to_cmap(sVals(sIdx), sClim, cmapPressureYellowRed(256));
        if useScatterNodes
            S = dataR_to_points2(ax, Rs, Xc(sIdx), Yc(sIdx));
            h.nodes_fill_S = scatter(ax, Xc(sIdx), Yc(sIdx), S, [1 1 1], 'filled','Marker','o','MarkerEdgeColor','none','HitTest','off');
            h.nodes_edge_S = scatter(ax, Xc(sIdx), Yc(sIdx), S*1.10, sRGB, 'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',0.75,'HitTest','off');
        else
            [h.nodes_fill_S, h.nodes_edge_S] = draw_node_circles_batched(ax, Xc(sIdx), Yc(sIdx), Rs, sRGB, [1 1 1]);
        end
    else
        h.nodes_fill_S = gobjects(1); h.nodes_edge_S = gobjects(1); sClim = [];
    end

    if ~isempty(dIdx)
        dClim = local_clim(dVals(dIdx));
        dRGB  = map_to_cmap(dVals(dIdx), dClim, cmapFlowTurquoiseBlue(256));
        if useScatterNodes
            S = dataR_to_points2(ax, Rd, Xc(dIdx), Yc(dIdx));
            h.nodes_fill_D = scatter(ax, Xc(dIdx), Yc(dIdx), S, [1 1 1], 'filled','Marker','o','MarkerEdgeColor','none','HitTest','off');
            h.nodes_edge_D = scatter(ax, Xc(dIdx), Yc(dIdx), S*1.10, dRGB, 'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',0.75,'HitTest','off');
        else
            [h.nodes_fill_D, h.nodes_edge_D] = draw_node_circles_batched(ax, Xc(dIdx), Yc(dIdx), Rd, dRGB, [1 1 1]);
        end
    else
        h.nodes_fill_D = gobjects(1); h.nodes_edge_D = gobjects(1); dClim = [];
    end

    setappdata(ax,'sClim',sClim);
    setappdata(ax,'dClim',dClim);

elseif strcmp(mode,'price') && hasX
    valsP = price_compact(:);
    p95 = prctile_safe(valsP,[5 95]);
    climPrice = [p95(1) p95(2)];
    if ~(isfinite(climPrice(1)) && isfinite(climPrice(2))) || climPrice(2) <= climPrice(1)
        lo = min(valsP); hi = max(valsP);
        if ~isfinite(lo) || ~isfinite(hi) || hi<=lo, lo = 0; hi = 1; end
        climPrice = [lo hi+eps];
    end
    nodeEdgeColors = map_to_cmap(valsP, climPrice, cmapPressureYellowRed(256));
    if useScatterNodes
        S = dataR_to_points2(ax, nodeR, Xc, Yc);
        h.nodes_fill = scatter(ax, Xc, Yc, S, [1 1 1], 'filled', 'Marker','o', 'MarkerEdgeColor','none','HitTest','off');
        h.nodes_edge = scatter(ax, Xc, Yc, S*1.10, nodeEdgeColors, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor','flat','LineWidth',0.75,'HitTest','off');
    else
        [h.nodes_fill, h.nodes_edge] = draw_node_circles_batched(ax, Xc, Yc, nodeR, nodeEdgeColors, [1 1 1]);
    end
    setappdata(ax,'climPrice',climPrice);

else
    % DEFAULT / SOLUTION / TOPOLOGY
    % Goal: non-gnodes -> PRESSURE colors; g-nodes -> FLOW colors
    applyPressureColors = hasX && ~isempty(pressure) && numel(pressure)==N && ~strcmp(mode,'flow') && ~isPreview;

    % 1) Base edge colors from PRESSURE for all nodes
    if applyPressureColors
        prP = prctile_safe(pressure(:),[5 95]);
        climP = [prP(1) prP(2)];
        if ~(isfinite(climP(1)) && isfinite(climP(2))) || climP(2)<=climP(1)
            climP = [min(pressure) max(pressure)+eps];
        end
        nodeEdgeColorsP = map_to_cmap(pressure(:),climP,cmapP);
    else
        nodeEdgeColorsP = repmat(colEdgeNode_default,N,1);
        if ~isempty(Gi), nodeEdgeColorsP(Gi,:) = repmat(colEdgeG, numel(Gi), 1); end
        climP = [];
    end

    % 2) Override ONLY g-nodes with FLOW colors (abs net at physical node)
    nodeEdgeColors = nodeEdgeColorsP;
    if doFlow && ~isempty(Gi) && ~isempty(climF)
        gvals = gnode_net_at_nodes(Gi);          % | q + G*xd - G*xs |
        colG  = map_to_cmap(gvals, climF, cmapF);% same cmap/limits as pipes/comps
        nodeEdgeColors(Gi,:) = colG;             % g-nodes show FLOW, others keep PRESSURE
    end

    % 3) Render: white fill, colored border
    if useScatterNodes
        S = dataR_to_points2(ax, nodeR, Xc, Yc);
        h.nodes_fill = scatter(ax, Xc, Yc, S, [1 1 1], 'filled', ...
            'Marker','o','MarkerEdgeColor','none','HitTest','off');
        h.nodes_edge = scatter(ax, Xc, Yc, S*1.10, nodeEdgeColors, 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor','flat', ...
            'LineWidth',0.75,'HitTest','off');
    else
        [h.nodes_fill, h.nodes_edge] = draw_node_circles_batched( ...
            ax, Xc, Yc, nodeR, nodeEdgeColors, [1 1 1]);
    end

    % Keep pressure clim for the pressure colorbar in "solution"
    setappdata(ax,'climP',climP);
end


%% ---------- draw compressor triangles ----------
if nC>0
    angles = atan2(Yc(ct)-Yc(cf), Xc(ct)-Xc(cf));
    if strcmp(mode,'pressure') || strcmp(mode,'gnode_flow') || strcmp(mode,'price')
        triRGB = repmat([0.8 0.8 0.8], nC, 1);
    elseif isPreview
        triRGB = repmat([0.67 0.13 1.0], nC, 1);
    elseif hasX && doFlow && ~isempty(flowC) && numel(flowC)==nC
        triRGB = map_to_cmap(abs(flowC(:)), climF, cmapF);
    else
        triRGB = repmat([0.67 0.13 1.0], nC, 1);
    end
    h.comp_tris = draw_triangles_patch(ax, triCenters, triSide, angles, triRGB);
else
    h.comp_tris = gobjects(0,1);
end

%% ---------- colorbars ----------
% Flow colorbar in 'solution'/'flow'
if (strcmp(mode,'flow') || strcmp(mode,'solution')) && doFlow && ~isPreview
    axpos  = get(ax,'Position'); gutter = 0.01; cbw = 0.02;
    axF = axes('Units','normalized', 'Position', axpos, 'Visible','off', 'Color','none');
    colormap(axF, cmapF);
    caxis(axF,climF);
    h.cbFlow = colorbar(axF, 'Location','eastoutside'); h.cbFlow.Units = 'normalized';
    pos = h.cbFlow.Position;
    pos(1)=axpos(1)+axpos(3)+gutter; pos(3)=cbw; pos(2)=axpos(2)+0.02*axpos(4); pos(4)=0.96*axpos(4);
    h.cbFlow.Position = pos; set(h.cbFlow,'FontSize',FS.cbTicks);
    ylabel(h.cbFlow,'$|$ Mass Flow $|$ (kg/s)','Interpreter','latex','fontsize',FS.cbLabel);
else
    h.cbFlow = gobjects(1);
end

% Pressure CB in 'solution'/'pressure'
climP = getappdata(ax,'climP');
showPressureCB = ~isempty(climP) && ~isPreview && ~strcmp(mode,'flow') && ~strcmp(mode,'gnode_flow') && ~strcmp(mode,'price');
if showPressureCB
    axpos  = get(ax,'Position'); gutter = 0.01; cbh = 0.03;
    axP = axes('Units','normalized', 'Position', axpos, 'Visible','off', 'Color','none');
    colormap(axP, cmapP); caxis(axP, climP);
    h.cbPress = colorbar(axP, 'Location','southoutside'); h.cbPress.Units = 'normalized';
    pos = h.cbPress.Position;
    pos(1)=axpos(1)+0.02*axpos(3); pos(3)=0.96*axpos(3); pos(2)=axpos(2)-gutter-cbh; pos(4)=cbh;
    h.cbPress.Position = pos; set(h.cbPress,'FontSize',FS.cbTicks);
    ticks = get(h.cbPress,'Ticks'); h.cbPress.Ruler.Exponent = 0;
    h.cbPress.TickLabels = arrayfun(@(t) sprintf('%.2f', t/1e6), ticks, 'UniformOutput', false);
    xlabel(h.cbPress,'Pressure (MPa)','Interpreter','latex','fontsize',FS.cbLabel);
else
    h.cbPress = gobjects(1);
end

% GNODE colorbars (flow):
if strcmp(mode,'gnode_flow') && hasX
    sClim = getappdata(ax,'sClim');  dClim = getappdata(ax,'dClim');

    if ~isempty(sClim)
        axpos = get(ax,'Position'); gutter=0.01; cbh=0.03;
        axS = axes('Units','normalized','Position',axpos,'Visible','off','Color','none');
        colormap(axS, cmapPressureYellowRed(256)); caxis(axS, sClim);
        h.cbSupply = colorbar(axS, 'Location','southoutside'); h.cbSupply.Units='normalized';
        pos = h.cbSupply.Position;
        pos(1)=axpos(1)+0.02*axpos(3); pos(3)=0.96*axpos(3); pos(2)=axpos(2)-gutter-cbh; pos(4)=cbh;
        h.cbSupply.Position = pos; set(h.cbSupply,'FontSize',FS.cbTicks);
        xlabel(h.cbSupply, '$|$ Supply Mass Flow $|$ (kg/s)', 'Interpreter','latex','fontsize',FS.cbLabel);
    else
        h.cbSupply = gobjects(1);
    end

    if ~isempty(dClim)
        axpos = get(ax,'Position'); gutter=0.01; cbw=0.02;
        axD = axes('Units','normalized','Position',axpos,'Visible','off','Color','none');
        colormap(axD, cmapFlowTurquoiseBlue(256)); caxis(axD, dClim);
        h.cbDemand = colorbar(axD, 'Location','eastoutside'); h.cbDemand.Units='normalized';
        pos = h.cbDemand.Position;
        pos(1)=axpos(1)+axpos(3)+gutter; pos(3)=cbw; pos(2)=axpos(2)+0.02*axpos(4); pos(4)=0.96*axpos(4);
        h.cbDemand.Position = pos; set(h.cbDemand,'FontSize',FS.cbTicks);
        ylabel(h.cbDemand, '$|$ Demand Mass Flow $|$ (kg/s)', 'Interpreter','latex','fontsize',FS.cbLabel);
    else
        h.cbDemand = gobjects(1);
    end
end

% PRICE colorbar:
if strcmp(mode,'price') && hasX
    climPrice = getappdata(ax,'climPrice');
    if ~isempty(climPrice)
        axpos  = get(ax,'Position'); gutter = 0.01; cbh = 0.03;
        axPrice = axes('Units','normalized','Position',axpos,'Visible','off','Color','none');
        colormap(axPrice, cmapPressureYellowRed(256)); caxis(axPrice, climPrice);
        h.cbPrice = colorbar(axPrice, 'Location','southoutside'); h.cbPrice.Units='normalized';
        pos = h.cbPrice.Position;
        pos(1)=axpos(1)+0.02*axpos(3); pos(3)=0.96*axpos(3); pos(2)=axpos(2)-gutter-cbh; pos(4)=cbh;
        h.cbPrice.Position = pos; set(h.cbPrice,'FontSize',FS.cbTicks);
        xlabel(h.cbPrice, 'Locational Marginal Price ($ \$ $/kg/s)', 'Interpreter','latex','fontsize',FS.cbLabel);
    else
        h.cbPrice = gobjects(1);
    end
end

%% ---------- frame ----------
[xmin, xmax, ymin, ymax] = compute_plot_bbox( ...
    Xc, Yc, nodeR, [], [], [], triCenters, triSide);
if ~isfinite(xmin) || ~isfinite(xmax) || xmin==xmax, xmin = min(Xc); xmax = max(Xc); end
if ~isfinite(ymin) || ~isfinite(ymax) || ymin==ymax, ymin = min(Yc); ymax = max(Yc); end
dx = xmax - xmin; dy = ymax - ymin;
if dx <= 0, dx = max(eps, range(Xc)); end
if dy <= 0, dy = max(eps, range(Yc)); end
pad = 0.01;
xlim(ax, [xmin - pad*dx, xmax + pad*dx]);
ylim(ax, [ymin - pad*dy, ymax + pad*dy]);

if strcmp(mode,'gnode_flow') || strcmp(mode,'price') || N > BIG_N
    set(findall(ax,'Type','patch'),'Clipping','on');
    set(findall(ax,'Type','line'), 'Clipping','on');
else
    set(findall(ax,'Type','patch'),'Clipping','off');
    set(findall(ax,'Type','line'), 'Clipping','off');
end

set(gcf,'PaperPositionMode','auto');
set(ax, 'FontSize', FS.axes);
if ~isempty(opt.ModelName)
    title(ax, opt.ModelName, 'Interpreter','latex', 'FontSize', FS.title);
end
drawnow limitrate nocallbacks;
hold(ax,'off');

%% ---------- local helper (scoping) ----------
    function vext = pad_or_pass_ext(v, maxid_, ext_ids_)
        v = v(:);
        if numel(v)==maxid_
            vext = v;
        elseif numel(v)==numel(ext_ids_)
            vext = zeros(maxid_,1); vext(ext_ids_) = v;
        else
            vext = zeros(maxid_,1);
            if numel(v)<=numel(ext_ids_)
                n = min(numel(v), numel(ext_ids_));
                vext(ext_ids_(1:n)) = v(1:n);
            end
        end
    end

end
% ===================== end render_one =====================


%% ===================== helpers =====================

function [hs, bins] = draw_grouped_segments(ax, X1, Y1, X2, Y2, C, w, nbins)
E = numel(X1);
if E==0
    hs = gobjects(0,1); bins = [];
    return;
end
edgesW = linspace(min(w), max(w)+eps, max(2,nbins)+1);
[~,bin] = histc(w, edgesW);
nb = max(bin); if isempty(nb) || nb<1, nb=1; end
hs = gobjects(nb,1);
for b = 1:nb
    sel = (bin==b);
    if any(sel)
        wbin = mean(edgesW([b b+1]));
        hs(b) = draw_segments_patch(ax, X1(sel),Y1(sel),X2(sel),Y2(sel), C(sel,:), wbin);
    else
        hs(b) = gobjects(1);
    end
end
bins = bin;
end

function hp = draw_segments_patch(ax, X1, Y1, X2, Y2, C, lw)
E = numel(X1);
if E==0
    hp = gobjects(1); return;
end
X1 = single(X1(:)); Y1 = single(Y1(:)); X2 = single(X2(:)); Y2 = single(Y2(:));
X = [X1.'; X2.']; Y = [Y1.'; Y2.']; Z = zeros(2,E,'single');
Cd = zeros(2, E, 3, 'single');
C = single(max(0,min(1,C)));
Cd(1,:,1) = C(:,1).';  Cd(1,:,2) = C(:,2).';  Cd(1,:,3) = C(:,3).';
Cd(2,:,:) = Cd(1,:,:);
hp = surface(ax, X, Y, Z, Cd, ...
    'FaceColor','none', 'EdgeColor','flat', 'CDataMapping','direct', ...
    'MeshStyle','column', 'LineWidth', lw, 'Marker','none', 'HitTest','off');
end

function hs = draw_segments_binned(ax, X1, Y1, X2, Y2, cVals, wVals, nColorBins, nWidthBins, cmap)
if isempty(X1), hs = gobjects(0); return; end
cVals = double(cVals(:));
wVals = double(wVals(:));
cLo = min(cVals); cRa = max(eps, range(cVals));
cIdx = 1 + floor( min(1, max(0,(cVals - cLo) / cRa)) * (nColorBins-1) );
wEdges = linspace(min(wVals), max(wVals)+eps, max(2,nWidthBins));
[~, wIdx] = histc(wVals, wEdges); wIdx(wIdx==0)=1;
hs = gobjects(nColorBins*nWidthBins,1); k = 0;

for ic = 1:nColorBins
    cRow = round( 1 + (ic-1)*(size(cmap,1)-1)/(nColorBins-1) );
    col  = cmap(cRow,:);
    for iw = 1:nWidthBins
        sel = (cIdx==ic) & (wIdx==iw);
        if any(sel)
            [XS, YS] = seg_nan([X1(sel) X2(sel)], [Y1(sel) Y2(sel)]);
            k = k+1;
            hs(k) = plot(ax, XS, YS, '-', ...
                'Color', col, ...
                'LineWidth', mean(wEdges([iw iw+1])), ...
                'HitTest','off');
        end
    end
end
hs = hs(1:k);
end

function [hfill, hedge] = draw_node_circles_batched(ax, Xc, Yc, R, edgeRGB, faceRGB)
m   = 32; t   = linspace(0,2*pi,m+1).'; t(end) = [];
N   = numel(Xc); ux  = cos(t); uy = sin(t);
V   = zeros(N*m,2,'single'); F   = reshape(1:(N*m), m, N).';
Xc = single(Xc); Yc = single(Yc); R = single(R);
for i=1:N
    idx = (i-1)*m + (1:m);
    V(idx,1) = Xc(i) + R(i)*ux;
    V(idx,2) = Yc(i) + R(i)*uy;
end
faceC_white = repmat(faceRGB, N, 1);
FV_white = zeros(N*m,3,'single');
for i=1:N
    idx = (i-1)*m + (1:m);
    FV_white(idx,:) = repmat(faceC_white(i,:), m, 1);
end
hfill = patch(ax, 'Faces',F, 'Vertices',V, 'FaceColor','flat', 'EdgeColor','none', ...
    'FaceVertexCData', FV_white, 'HitTest','off');
FV_edge = zeros(N*m,3,'single');
edgeRGB = single(edgeRGB);
for i=1:N
    idx = (i-1)*m + (1:m);
    FV_edge(idx,:) = repmat(edgeRGB(i,:), m, 1);
end
hedge = patch(ax, 'Faces',F, 'Vertices',V, ...
    'FaceColor','none', 'EdgeColor','flat', 'LineWidth',1.0, ...
    'FaceVertexCData', FV_edge, 'HitTest','off');
end

function S = dataR_to_points2(ax, R_data, Xc, Yc)
% Convert data-units radius to scatter SizeData (points^2)
drawnow;  % ensure positions are realized
pp = getpixelposition(ax);
xlim_ = xlim(ax); ylim_ = ylim(ax);

px_per_x = pp(3) / max(eps, diff(xlim_));
px_per_y = pp(4) / max(eps, diff(ylim_));
px_per_data = min(px_per_x, px_per_y);
pt_per_px = 72/get(0,'ScreenPixelsPerInch');
pt_per_data = px_per_data * pt_per_px;

S = (2 * R_data * pt_per_data).^2;

% --- Sanitize for scatter ---
S(~isfinite(S)) = NaN; 
S(S<=0) = 1;
S = double(S);
end

function htri = draw_triangles_patch(ax, centers, sides, angles, edgeRGB)
M  = size(centers,1);
if M==0
    htri = gobjects(0,1); return;
end
m3 = 3; V  = zeros(M*m3,2,'single'); F  = zeros(M,m3);
centers = single(centers); sides = single(sides); angles = single(angles);
for k=1:M
    s   = sides(k);
    h   = s*sqrt(3)/2;
    loc = [  2*h/3,  -h/3,   -h/3;
        0   ,  -s/2,    s/2 ];
    R2  = single([cos(angles(k)) -sin(angles(k)); sin(angles(k)) cos(angles(k))]);
    tri = (R2*loc) + centers(k,:).';
    idx = (k-1)*m3 + (1:m3);
    V(idx,:) = tri.'; F(k,:) = idx;
end
FV_edge = zeros(M*m3,3,'single');
edgeRGB = single(edgeRGB);
for k=1:M
    idx = (k-1)*m3 + (1:m3);
    FV_edge(idx,:) = repmat(edgeRGB(k,:), m3, 1);
end
htri = patch(ax,'Faces',F,'Vertices',V, 'FaceColor',[1 1 1], 'EdgeColor','flat', ...
    'LineWidth',1.0, 'FaceVertexCData', FV_edge, 'HitTest','off');
end

function nodeR = node_radii_knn_fast(X, Y, k)
if nargin<3, k=1; end
XYs = single([X(:), Y(:)]);
N = size(XYs,1);
try
    Mdl = KDTreeSearcher(XYs,'Distance','euclidean');
    [~, D] = knnsearch(Mdl, XYs, 'K', k+1);
    nn = double(D(:, 1+k));
    safety=0.98; r_noOL=0.5*safety*nn;
catch
    rx=max(X)-min(X); ry=max(Y)-min(Y); boxA=max(rx*ry,eps);
    r_noOL = 0.5*sqrt((0.06*boxA)/(N*pi))*ones(N,1);
end
rx=max(X)-min(X); ry=max(Y)-min(Y); boxA=max(rx*ry,eps);
r_packcap=0.5*sqrt((0.06*boxA)/(N*pi))*ones(N,1);
dataSpan=max([rx,ry,eps]); dataMin=0.006*dataSpan; dataMax=0.10*dataSpan;
nodeR=min([r_noOL,r_packcap,dataMax*ones(N,1)],[],2);
nodeR=max(nodeR, min(dataMin, r_noOL));
nodeR=min(1.20*nodeR, r_noOL);

% --- Tiny positive lower bound in data units ---
span = max([max(X)-min(X), max(Y)-min(Y), eps]);
floorR = 1e-6 * span;           % very small, relative to plot size
nodeR(~isfinite(nodeR)) = floorR;
nodeR = max(nodeR, floorR);
end

function [centers, side, triMinSide] = compressor_triangles_adaptive_fast( ...
    X, Y, cf, ct, nodeR, Gi, maxIter, tri_scale, tri_q, simpleMode)
M = numel(cf);
centers = zeros(M,2);
side    = zeros(M,1);
triMinSide = NaN;
if M==0, return; end

mx = 0.5*(X(cf) + X(ct));
my = 0.5*(Y(cf) + Y(ct));
centers = [mx my];
L  = hypot(X(ct) - X(cf), Y(ct) - Y(cf));

if simpleMode
    s_len = 0.60 * L;
    if ~isempty(Gi)
        rG = nodeR(Gi); rG = rG(rG>0);
        if isempty(rG), rG = nodeR(nodeR>0); end
    else
        rG = nodeR(nodeR>0);
    end
    rGmean = mean(rG);
    if ~isfinite(rGmean) || isempty(rGmean), rGmean = 0.01; end
    s_target  = 2.0 * rGmean;
    s_floor   = 0.9 * s_target;
    side = tri_scale * max(s_len, s_floor);
    % --- enforce max allowed compressor size in simple mode
    s_cap = tri_scale * sqrt(3) * max(nodeR);   % same cap logic as detailed path
    side  = min(side, s_cap);

    triMinSide = s_floor;
    return;
end

R_other = zeros(M,1);
XYs = single([X(:),Y(:)]);
try
    Mdl = KDTreeSearcher(XYs,'Distance','euclidean');
    kNN = 16;
    for m = 1:M
        mid = single([mx(m) my(m)]);
        [idx, D] = knnsearch(Mdl, mid, 'K', kNN);
        idx = idx(:);
        idx(idx==cf(m) | idx==ct(m)) = [];
        if isempty(idx)
            R_other(m) = inf;
        else
            d_to_nodes = double(D(1:numel(idx)))' - nodeR(idx);
            R_other(m) = max(0, 0.98 * min(d_to_nodes));
        end
    end
catch
    R_other(:) = inf;
end
r12 = max(nodeR(cf), nodeR(ct));
R_ends = max(0, 0.98*(0.5*L - r12));
R_clear = min(R_other, R_ends);
s_from_clear = sqrt(3) * R_clear;
s_from_len = 0.95 * L;

if ~isempty(Gi)
    rG = nodeR(Gi); rG = rG(rG>0);
    if isempty(rG), rG = nodeR(nodeR>0); end
else
    rG = nodeR(nodeR>0);
end
rGmean = mean(rG);
if ~isfinite(rGmean) || isempty(rGmean), rGmean = mean(nodeR(nodeR>0)); end
if ~isfinite(rGmean) || isempty(rGmean), rGmean = 0.01; end
s_target  = 2.0 * rGmean;
s_floor   = 0.9 * s_target;
triMinSide = s_floor;

s_hardcap = sqrt(3) * max(nodeR) * ones(M,1);
s0 = min([s_from_clear, s_from_len, s_hardcap], [], 2);
s0 = tri_scale * s0;
side = max(s0, min(s_floor, max(s_hardcap))*ones(M,1));

R = side ./ sqrt(3);
XY = centers;
try
    for it = 1:maxIter
        changed = false;
        Mdl2 = KDTreeSearcher(single(XY));
        idxs = rangesearch(Mdl2, single(XY), double(R)+double(R));
        for a = 1:M
            nb = idxs{a}; nb(nb<=a) = [];
            for b = nb
                d = hypot(XY(a,1)-XY(b,1), XY(a,2)-XY(b,2));
                if R(a) + R(b) > d && d>0
                    excess = R(a) + R(b) - d;
                    take   = 0.5 * excess;
                    Ra = max(R(a) - take, s_floor/sqrt(3));
                    Rb = max(R(b) - take, s_floor/sqrt(3));
                    if (Ra < R(a)) || (Rb < R(b))
                        R(a) = Ra; R(b) = Rb; changed = true;
                    end
                end
            end
        end
        if ~changed, break; end
    end
catch
    % keep current R
end
side = sqrt(3) * R;
end

function [xmin, xmax, ymin, ymax] = compute_plot_bbox(Xc, Yc, nodeR, slack_x, slack_y, slack_r, triCenters, triSide)
xmin = min(Xc - nodeR);  xmax = max(Xc + nodeR);
ymin = min(Yc - nodeR);  ymax = max(Yc + nodeR);
if ~isempty(slack_r)
    xmin = min(xmin, min(slack_x - slack_r));
    xmax = max(xmax, max(slack_x + slack_r));
    ymin = min(ymin, min(slack_y - slack_r));
    ymax = max(ymax, max(slack_y + slack_r));
end
if ~isempty(triSide)
    Rtri = triSide ./ sqrt(3);
    xmin = min(xmin, min(triCenters(:,1) - Rtri));
    xmax = max(xmax, max(triCenters(:,1) + Rtri));
    ymin = min(ymin, min(triCenters(:,2) - Rtri));
    ymax = max(ymax, max(triCenters(:,2) + Rtri));
end
end

function [XS,YS]=seg_nan(X2,Y2)
n=size(X2,1);
XS=reshape([X2.'; nan(1,n)],[],1);
YS=reshape([Y2.'; nan(1,n)],[],1);
end

function C = map_to_cmap(v,clim,cmap)
v=double(v); lo=clim(1); hi=clim(2);
if hi<=lo, hi=lo+eps; end
t=(v-lo)/(hi-lo); t=max(0,min(1,t));
idx=1+floor(t*(size(cmap,1)-1));
C=cmap(idx,:);
end

function p = prctile_safe(x,q)
x=sort(x(:)); if isempty(x), p=NaN(size(q)); return; end
q=max(0,min(100,q(:))); N=numel(x); p=zeros(size(q));
for ii=1:numel(q)
    if q(ii)==0, p(ii)=x(1); continue; end
    if q(ii)==100, p(ii)=x(end); continue; end
    pos=1+(N-1)*q(ii)/100; k=floor(pos); d=pos-k;
    if k>=N, p(ii)=x(end); else, p(ii)=x(k)*(1-d)+x(k+1)*d; end
end
end

function clim = local_clim(vals)
if isempty(vals), clim = []; return; end
p = prctile_safe(vals,[5 95]);
lo = 0; hi = max(p(2), max(vals)*1e-6);
clim = [lo hi];
end

function cmap = cmapFlowTurquoiseBlue(n)
if nargin < 1, n = 256; end
c1 = [0.6 1.0 1.0];
c2 = [0.0 0.0 0.4];
cmap = interp1([0 1], [c1; c2], linspace(0,1,n));
end

function cmap = cmapPressureYellowRed(n)
if nargin < 1, n = 256; end
c1 = [1.0 1.0 0.0];
c2 = [0.8 0.0 0.0];
cmap = interp1([0 1], [c1; c2], linspace(0,1,n));
end