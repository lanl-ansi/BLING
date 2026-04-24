function [net, bc, param] = read_input(folder)
% Load CSV inputs into net, bc, param
 
% Get list of all CSV files in the folder
files = dir(fullfile(folder, '*.csv'));

% Initialize network, boundary condition, and parameter structures
net = struct();
bc = struct();
param = struct();

% Read input data into array structures
array_structs();

% net.comps (after column removel): [comp label, from, to, min ratio, max ratio, min flow, max flow]
if isfield(net,'comps') && ~isempty(net.comps) && size(net.comps,2) >= 7
    net.comps(:,[2,7]) = [];
end

% net.pipes (after column removel): [pipe label, from, to, diameter, length, friction]
if isfield(net,'pipes') && ~isempty(net.pipes) && size(net.pipes,2) >= 2
    net.pipes(:,2) = [];
end

% net.nodes (after column removel): [node label, x, y, pmin, pmax, slackBool]
if isfield(net,'nodes') && ~isempty(net.nodes) && size(net.nodes,2) >= 8
    net.nodes(:,[2,7,8]) = [];
end

% net.gnodes (after column removel): [gnode label, physical node label]
if isfield(net,'gnodes') && ~isempty(net.gnodes) && size(net.gnodes,2) >= 2
    net.gnodes(:,2) = [];
end

% ---- Sort rows into columns
if isfield(net,'nodes')  && ~isempty(net.nodes),  net.nodes  = sortrows(net.nodes,  1); end
if isfield(net,'comps')  && ~isempty(net.comps),  net.comps  = sortrows(net.comps,  1); end
if isfield(net,'pipes')  && ~isempty(net.pipes),  net.pipes  = sortrows(net.pipes,  1); end
if isfield(net,'gnodes') && ~isempty(net.gnodes), net.gnodes = sortrows(net.gnodes, 1); end

if isfield(bc,'cd')   && ~isempty(bc.cd),   bc.cd   = sortrows(bc.cd',   1)'; end
if isfield(bc,'cs')   && ~isempty(bc.cs),   bc.cs   = sortrows(bc.cs',   1)'; end
if isfield(bc,'gbar') && ~isempty(bc.gbar), bc.gbar = sortrows(bc.gbar', 1)'; end
if isfield(bc,'qbar') && ~isempty(bc.qbar), bc.qbar = sortrows(bc.qbar', 1)'; end
if isfield(bc,'smax') && ~isempty(bc.smax), bc.smax = sortrows(bc.smax', 1)'; end
if isfield(bc,'dmax') && ~isempty(bc.dmax), bc.dmax = sortrows(bc.dmax', 1)'; end

% ---- Relabel IDs from 1 to N
if isfield(net,'nodes') && ~isempty(net.nodes)
    net.nodes(:,1) = (1:size(net.nodes,1)).';
end
if isfield(net,'comps') && ~isempty(net.comps)
    net.comps(:,1) = (1:size(net.comps,1)).';
end
if isfield(net,'pipes') && ~isempty(net.pipes)
    net.pipes(:,1) = (1:size(net.pipes,1)).';
end

% ================== read all CSVs ==================
    function array_structs()
        for i = 1:numel(files)
            fileName = files(i).name;
            filePath = fullfile(folder, fileName);
            baseName = matlab.lang.makeValidName(erase(fileName, '.csv'));

            % Safely read first two lines
            [firstLine, secondLine] = first_two_lines(filePath);

            % If file is empty or unreadable, skip
            if isempty(firstLine) || isempty(secondLine)
                warning('read_input:EmptyOrUnreadable','Skipping "%s" (empty or unreadable).', fileName);
                continue
            end

            % ---- Detect file type
            secondRowCells = strsplit(secondLine, ',');
            isParam = numel(secondRowCells) == 2 ...
                      && isnan(str2double(secondRowCells{1})) ...
                      && ~isnan(str2double(secondRowCells{2}));

            if isParam
                % param: (text header in col 1, number in col 2)
                tbl = readtable(filePath, 'ReadVariableNames', false);
                if size(tbl,2) < 2
                    warning('read_input:ParamFormat','Param file "%s" has <2 columns. Skipping.', fileName);
                    continue
                end
                for k = 1:height(tbl)
                    % parameter type
                    key = tbl{k,1};
                    if iscell(key), key = key{1}; end
                    key = string(key);
                    key = matlab.lang.makeValidName(key);
                    % parameter value
                    value = tbl{k,2};
                    if iscell(value), value = value{1}; end
                    if ischar(value) || isstring(value)
                        value = str2double(string(value));
                    end
                    if ~isscalar(value) || ~isnumeric(value)
                        warning('read_input:ParamValue','Param "%s" in "%s" not scalar numeric. Skipping row.', key, fileName);
                        continue
                    end
                    param.(char(key)) = value;
                end

            else
                % Determine if first row is numeric -> bc (no headers)
                firstRowCells = strsplit(firstLine, ',');
                isNumericRow = all(cellfun(@(x) ~isnan(str2double(x)), firstRowCells));

                if isNumericRow
                    % bc: no headers, all numeric
                    matrixData = readmatrix(filePath);
                    if ~isnumeric(matrixData)
                        warning('read_input:BCNonNumeric','BC file "%s" is not numeric. Skipping.', fileName);
                        continue
                    end
                    bc.(baseName) = matrixData;
                else
                    % net: has headers in first row
                    matrixData = readmatrix(filePath, 'NumHeaderLines', 1);
                    if ~isnumeric(matrixData)
                        warning('read_input:NetNonNumeric','Net file "%s" is not numeric. Skipping.', fileName);
                        continue
                    end
                    net.(baseName) = matrixData;
                end
            end
        end
    end

% ---- Helper
    function [firstLine, secondLine] = first_two_lines(filePath)
        firstLine = ''; secondLine = '';
        fid = fopen(filePath,'r');
        if fid < 0, return; end
        c = onCleanup(@() fclose(fid));
        firstLine  = fgetl(fid); if ~ischar(firstLine),  firstLine  = ''; return; end
        secondLine = fgetl(fid); if ~ischar(secondLine), secondLine = ''; return; end
    end
end