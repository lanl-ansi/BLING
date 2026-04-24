function output(x, par, net, param, output_folder, history)

% Define data names and indices for the CSV files
if param.save_data == 1
    names = {'pipe_flow_in', ... % size: |pipes| -- units: kg/s
        'pipe_flow_out', ...     % size: |pipes| -- units: kg/s
        'comp_pressure_in', ...  % size: |comps| -- units: Pa
        'comp_pressure_out', ... % size: |comps| -- units: Pa
        'comp_flow_in', ...      % size: |comps| -- units: kg/s
        'comp_flow_out', ...     % size: |comps| -- units: kg/s
        'comp_power', ...        % size: |comps| -- units: J/s (W) or MJ/s 
        'node_pressures', ...    % size: |nodes| -- units: Pa
        'slack_flows', ...       % size: |slack| -- units: kg/s
        'gNode_inflows', ...     % size: |gnode| -- units: kg/s
        'gNode_outflows', ...    % size: |gnode| -- units: kg/s
        'comp_ratios', ...       % size: |comps| -- units: n/a
        'LMP_gnode', ...         % size: |gnode| -- units: Dollars($) / (kg/s)
        'LMP_nodes', ...         % size: |nodes| -- units: Dollars($) / (kg/s)
        'nodal_flow'};           % size: |nodes| -- units: kg/s

    sets = {'pipes', ...
        'pipes', ...
        'nodes', ...
        'nodes', ...
        'comps', ...
        'comps', ...
        'nodes', ...
        'nodes', ...
        'comps', ...
        'gnode_s', ...
        'gnode_d', ...
        'nodes', ...
        'gnode_s', ...
        'nodes', ...
        'nodes'};

    if exist('par.Cp', 'var') == 0
        par.Cp = (par.C > 0);
        par.Cm = (par.C < 0);
    end

    for i = 1:numel(names)
        % Compute subset data
        index = par.index.(sets{i});
        data = compute_data(names{i}, index);

        % Define the filename using the predefined name
        filename = fullfile(output_folder, [names{i} '.csv']);

        % Write the subset data to a CSV file
        writematrix(data, filename);
    end
end

% --- Ouptput sparsity pattern
if param.sparse_pattern == 1
    write_sparsity_csv(par, net);
end

    function data = compute_data(name, index)
        % Explicit data
        data = [index'-min(index)+1; x(index)'];

        % Implicit data
        if name == "comp_pressure_in"
            dat = par.Cm * x(index);
            data = [1:par.n.comps; dat'];
        end
        if name == "comp_pressure_out"
            dat = par.Cp * x(index);
            data = [1:par.n.comps; dat'];
        end
        if name == "comp_power"
            r = param.specific_heat_capacity_ratio;
            ratio = (r - 1) / r;
            a = 286.76 * param.Temperature / (ratio * param.Gas_specific_gravity); % Adiabatic work coefficient
            dat = a * x(par.index.comps) .* max( ( (par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes)) ).^ratio - 1,0);
            data = [1:par.n.comps; dat' / 10^6]; % Units: MW
        end
        if name == "slack_flows"
            mask = (sum(par.Cm * par.S', 2)) ~= 0;
            dat = x(index(mask));

            n = numel(dat);
            data = [1:n; dat(:).'];
        end
        if name == "comp_ratios"
            dat = (par.Cp * x(index)) ./ (par.Cm * x(index));
            index = par.index.comps;
            data = [index'-min(index)+1; dat'];
        end
        if name == "LMP_gnode"
            dat = history(end).lambda_gnode;
            data = [index'-min(index)+1; dat'];
        end
        if name == "LMP_nodes"
            dat = history(end).lambda_all;
            data = [index'-min(index)+1; dat'];
        end
        if name == "nodal_flow"
            index = [par.index.comps; par.index.pipes];
            dat =  [par.C', par.P'] * x(index);
            index = par.index.nodes;
            data = [index'-min(index)+1; dat'];
        end
    end

    function write_sparsity_csv(par, net)
        % Slack pressure rows vs nodes: Abase(row, par.index.nodes) = par.S
        row_slack = (1:par.n.slack).';
        col_nodes = par.index.nodes(:);
        [iS, jS]  = find(par.S);
        rows_S = row_slack(iS);
        cols_S = col_nodes(jS);

        % Nodal balance block rows:
        row_nbal = par.n.slack + par.n.pipes + (1:par.n.nodes).';

        % Abase(row_nbal, par.index.comps) = par.C.'
        [iC, jC] = find(par.C.');
        rows_C = row_nbal(iC);
        cols_C = par.index.comps(jC);

        % Abase(row_nbal, par.index.pipes) = par.P.'
        [iP, jP] = find(par.P.');
        rows_P = row_nbal(iP);
        cols_P = par.index.pipes(jP);

        % Abase(row_nbal, par.index.gnode_s) = par.G
        [iG, jG] = find(par.G);
        rows_Gs = row_nbal(iG);
        cols_Gs = par.index.gnode_s(jG);

        % Abase(row_nbal, par.index.gnode_d) = -par.G
        rows_Gd = row_nbal(iG);
        cols_Gd = par.index.gnode_d(jG);

        fixed_sparsity = [rows_S, cols_S;
            rows_C, cols_C;
            rows_P, cols_P;
            rows_Gs, cols_Gs;
            rows_Gd, cols_Gd];

        % Deduplicate
        if ~isempty(fixed_sparsity)
            fixed_sparsity = unique(fixed_sparsity, 'rows');
        end

        % Reconstruct exactly the index pairs used in nonlinear_jacobian:
        pipeRows = par.n.slack + (1:par.n.pipes).';
        pipeCols = par.index.pipes(:);

        % PipeRow / PipeCol for node-pressure nonlinear terms
        PipeRow = par.n.slack + [net.pipes(:,1); net.pipes(:,1)];
        PipeCol_raw = [net.pipes(:,2); net.pipes(:,3)];

        PipeCol = PipeCol_raw;

        varying_sparsity = [PipeRow, PipeCol;
            pipeRows, pipeCols];

        % Deduplicate
        if ~isempty(varying_sparsity)
            varying_sparsity = unique(varying_sparsity, 'rows');
        end

        all_sparsity = unique([fixed_sparsity; varying_sparsity], 'rows');

        % Write CSVs
        prefix = 'parA';
        filename = fullfile(output_folder, [prefix '_fixed.csv']);
        writematrix(fixed_sparsity,   filename);
        filename = fullfile(output_folder, [prefix '_varying.csv']);
        writematrix(varying_sparsity, filename);
        filename = fullfile(output_folder, [prefix '_all.csv']);
        writematrix(all_sparsity,     filename);
    end
end