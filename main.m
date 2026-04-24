close all; clc

% --- Choose the network model to optimize
% model = 'LANL-8';
 model = 'LANL-30';

% --- Define input and output folder paths
input_folder = fullfile(pwd, 'examples', model, 'Input');
output_folder = fullfile(pwd, 'examples', model, 'Output');

% --- Read data from CSV files and create the structures net, bc, and param
[net, bc, param] = read_input(input_folder);

% --- Define objective and constraint matrices
[par, x0, history] = parameters(net, bc, param);

% --- Find optimal csv params
if param.suggest_params == 1
param = tune_csv_params(net, par, param, x0);
end

% --- SLP solve
tic
[x, history, par] = solveSLP(net, par, param, x0, history);

% --- Newton solve
if param.newton_solve == 1
    [x, par, history] = newton_solve(x, par, param, net, history, false);
end
history(end).comptime = toc;

% --- Rescale back to dimensional units
[x, x0, par] = rescaling(x, x0, par);

% --- Print solution properties

print_data(x, x0, par, history, model, net, param, input_folder, output_folder)


