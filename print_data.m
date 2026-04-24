function print_data(x, x0, par, history, model, net, param, input_folder, output_folder)

% Display solution properties on command window
fprintf('Elapsed time (sec): %.16f\n',history(end).comptime);
fprintf('Number of iterations for %s: %.0f\n', model, numel(history));
fprintf('Nodes: %.0f\n', par.n.nodes);
fprintf('Comps: %.0f\n', par.n.comps);
fprintf('Pipes: %.0f\n', par.n.pipes);
fprintf('Gnode: %.0f\n', par.n.gnodes);
fprintf('Slack: %.0f\n', par.n.slack);
fprintf('Length (mi): %.1f\n', 0.621371 * sum(net.pipes(:,5))/10^3);
fprintf('\n--- SOLUTION PROPERTIES ---\n');

fprintf('Residual = %.16f\n', history(end).constraint_violation);
fprintf('Upper pressure difference: MIN (p^max - p)./p * 100 = %.16f\n', 100 * min( (par.xmax(par.index.nodes) - x(par.index.nodes)) ./ x(par.index.nodes) ));
fprintf('Lower pressure difference: MIN (p - p^min)./p * 100 = %.16f\n', 100 * min( (x(par.index.nodes) - par.xmin(par.index.nodes)) ./ x(par.index.nodes) ));
fprintf('Up pressure viol nodes: # (p^max - p)./p<-.05 = %d\n', sum((par.xmax(par.index.nodes) - x(par.index.nodes))./ x(par.index.nodes)<-0.05)  );
fprintf('Low pressure viol nodes: # (p - p^min)./p<-.05 = %d\n', sum((x(par.index.nodes) - par.xmin(par.index.nodes))./ x(par.index.nodes)<-0.05)  );

fprintf('Max compress ratio: MAX u = %.16f\n', max( ( (par.Cp * x(par.index.nodes)) ./ (par.Cm *x( par.index.nodes)) ) ));
fprintf('Avg compress ratio: AVG u = %.16f\n', mean( ( (par.Cp * x(par.index.nodes)) ./ (par.Cm *x( par.index.nodes)) ) ));
fprintf('Min compress ratio: MIN u = %.16f\n', min( ( (par.Cp * x(par.index.nodes)) ./ (par.Cm *x( par.index.nodes)) ) ));

fprintf('\n--- SOLUTION NOMINATIONS ---\n');

fprintf('Supply ratio: SUM s/s^max = %.16f\n', sum( x(par.index.gnode_s)) / sum(par.xmax(par.index.gnode_s)));
fprintf('Demand ratio: SUM d/d^max = %.16f\n', sum(x(par.index.gnode_d) + par.G' * par.q)/sum(par.xmax(par.index.gnode_d) + par.G' * par.q));
fprintf('Total delivered flow kg/s = %.16f\n', sum(x(par.index.gnode_d) + par.G' * par.q));
fprintf('Total delivered energy MW = %.16f\n', sum(x(par.index.gnode_d) + par.G' * par.q)/0.01898);
fprintf('Total delivered energy mmscfd = %.16f\n', sum(x(par.index.gnode_d) + par.G' * par.q)/0.2404);
r = param.specific_heat_capacity_ratio;
ratio = (r - 1) / r;
a = 286.76 * param.Temperature / (ratio * param.Gas_specific_gravity); % Adiabatic work coefficient
energy = sum(a * x(par.index.comps) .* max(( (par.Cp * x(par.index.nodes)) ./ (par.Cm * x(par.index.nodes))).^ratio - 1,0));
fprintf('Comps energy: SUM engy MW = %.16f\n', energy /10^6); % units: MW
fprintf('Total fuel factor percent = %.16f\n', energy /10^6/(sum(x(par.index.gnode_d) + par.G' * par.q)/0.01898)*100);

% --- Save data in CSV files located in output_folder
if param.save_data == 1 || param.sparse_pattern == 1
    output(x, par, net, param, output_folder, history)
end

% --- Plot network topology and solution overlay
lambda = history(end).lambda_all;
if (param.pressure_flow_plot == 1 || param.pressure_plot == 1 || param.flow_plot == 1 || ...
        param.gnode_flow_plot == 1 || param.LMP_plot == 1 )

    plot_network(net, par, param, lambda, x, 'num_width_bins',1,'flow_width', [1.5 1.5], 'ModelName', model);
end
if param.topology_only_plot == 1
    plot_network(net, par, param, lambda, 'num_width_bins',1,'flow_width', [1.5 1.5], 'ModelName', model);
end

% --- Update param.csv
if param.suggest_params == 1
    writetable(rows2vars(struct2table(param)), fullfile(input_folder, 'param.csv'), ...
        'WriteVariableNames', false, 'WriteRowNames', false);
end

% --- Plot convergence of SLP
if param.plot_convergence == 1
    convergence_plot(history, x0, x, par);
end
end

