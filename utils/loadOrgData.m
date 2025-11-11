function [Data, Encoding, Lower, Upper] = LoadOrgData(config)
% LOADORGDATA - load and prepare original problem data and bounds
%
% [Data, Encoding, Lower, Upper] = LoadOrgData(config)
%
% Reads problem data according to config.Pro and returns:
%  Data     - problem-specific raw data (format depends on problem)
%  Encoding - vector describing encoding type per decision variable (1 == numeric)
%  Lower    - vector of lower bounds for decision variables
%  Upper    - vector of upper bounds for decision variables
%
% The function is structured so adding new problem types is simple:
%   - implement a local helper load_<PRONAME> that returns the four outputs
%   - Use `config.Pro = "NEWPRO"` to select it

    % Normalize the problem name to a char row for robust comparisons
    pro = upper(char(config.Pro));

    switch pro
        case 'TSP'
            [Data, Encoding, Lower, Upper] = load_TSP();
        case 'KP'
            [Data, Encoding, Lower, Upper] = load_KP();
        case 'MOTSP'
            [Data, Encoding, Lower, Upper] = load_MOTSP();
        case 'MOKP'
            [Data, Encoding, Lower, Upper] = load_MOKP();
        case 'PO'
            [Data, Encoding, Lower, Upper] = load_PO();
        case 'OPF'
            [Data, Encoding, Lower, Upper] = load_OPF();
        case 'CRT'
            [Data, Encoding, Lower, Upper] = load_CRT(config);
        otherwise
            error('LoadOrgData:UnknownProblem', 'please input Pro = TSP or KP or MOTSP or MOKP or PO or OPF or CRT');
    end

end

%% --- Per-problem loaders ---
function [Data, Encoding, Lower, Upper] = load_TSP()
    dat = load('Data_TSP-D200.mat');
    Data = dat.Data;
    Encoding = ones(1, numel(Data));
    Lower = -reshape(Data, 1, []);
    Upper = 1 - reshape(Data, 1, []);
end

function [Data, Encoding, Lower, Upper] = load_KP()
    % Knapsack: weights and profits provided in separate .mat files
    tmp = load('P-100D.mat','P'); P = tmp.P; clear tmp
    tmp = load('W-100D.mat','W'); W = tmp.W; clear tmp
    Data = [P, W];
    Encoding = ones(1, numel(Data));
    Lower = -reshape(Data, 1, []);
    Upper = 1 - reshape(Data, 1, []);
end

function [Data, Encoding, Lower, Upper] = load_MOTSP()
    tmp = load('MOTSP-D90.mat','C');
    Data = tmp.C;
    % Encoding size is application-specific; preserve original expression
    Encoding = ones(1, 2*10*(10-1)/2);
    index1 = find(tril(Data{1}, -1) ~= 0);
    index2 = find(tril(Data{2}, -1) ~= 0);
    temp1 = tril(Data{1}, -1);
    temp2 = tril(Data{2}, -1);
    Lower = -[temp1(index1); temp2(index2)]';
    Upper = 1 - [temp1(index1); temp2(index2)]';
end

function [Data, Encoding, Lower, Upper] = load_MOKP()
    tmp = load('MOP-M2-D250.mat','MOP'); MOP = tmp.MOP; clear tmp
    tmp = load('MOW-M2-D250.mat','MOW'); MOW = tmp.MOW; clear tmp
    Data = [MOP; MOW];
    Data = reshape(Data, 1, []);
    Encoding = ones(1, numel(Data));
    Lower = -reshape(Data, 1, []);
    Upper = 1 - reshape(Data, 1, []);
end

function [Data, Encoding, Lower, Upper] = load_PO()
    tmp = load('Dataset_PO.mat','Dataset'); Dataset = tmp.Dataset; clear tmp
    Data = Dataset.('data1000');
    Data = Data(1:20, 1:10);
    Datamin = min(Data, [], 'all');
    Datamax = max(Data, [], 'all');
    Datanorm = (Data - Datamin) / (Datamax - Datamin);
    Encoding = ones(1, numel(Datanorm));
    Lower = -reshape(Datanorm, 1, []);
    Upper = 1 - reshape(Datanorm, 1, []);
end

function [Data, Encoding, Lower, Upper] = load_OPF()
    tmp = load('mpc.mat');
    Data = tmp.mpc.gencost(:, 5:6);
    max_1 = max(Data(:,1)); min_1 = min(Data(:,1));
    max_2 = max(Data(:,2)); min_2 = min(Data(:,2));
    Datanorm(:,1) = (Data(:,1) - min_1) / (max_1 - min_1);
    Datanorm(:,2) = (Data(:,2) - min_2) / (max_2 - min_2);
    Encoding = ones(1, numel(Data));
    Lower = -reshape(Datanorm, 1, []);
    Upper = 1 - reshape(Datanorm, 1, []);
end

function [Data, Encoding, Lower, Upper] = load_CRT(config)
    % CRT: sample zero entries to reach a fixed size and return only
    % the non-zero subset (preserve original behavior)
    tmp = load('Eyeball_L_PoiIndex_weight.mat');
    Eyeball_L_PoiIndex_weight = tmp.Eyeball_L_PoiIndex_weight; clear tmp
    Eyeball_L = Eyeball_L_PoiIndex_weight;
    nonzp = find(Eyeball_L ~= 0);
    zp = find(Eyeball_L == 0);
    num_to_select = 500 - numel(nonzp);
    if num_to_select > 0
        szp = randsample(zp, num_to_select);
    else
        szp = [];
    end
    sp = [nonzp; szp]; %#ok<NASGU>
    Data = Eyeball_L(nonzp);
    Datanorm = (Data - min(Data)) / (max(Data) - min(Data));
    Encoding = ones(1, numel(Data));
    Lower = -config.scalefactor * Datanorm';
    Upper = config.scalefactor * (1 - Datanorm');
    % save the non-zero indices into the PlatEMO folder as previously
    try
        root = fileparts(fileparts(mfilename('fullpath')));
        save(fullfile(root, 'PlatEMO 4.0', 'sp.mat'), 'nonzp');
    catch
        % ignore save failures to avoid breaking data load
    end
end

