% ORIGIN - compute and cache the original (baseline) performance for a problem
%
% origin = Origin(Pro, Data, A_attacked, n, maxfe, run, save_foder, varargin)
%
% Inputs:
%   Pro         - problem identifier string: 'TSP','KP','MOTSP','MOKP','PO','OPF','CRT'
%   Data        - problem-specific data (format varies by problem)
%   A_attacked  - name of the attacked algorithm constructor (string). The
%                 function is constructed with str2func and called with varargin
%   n           - population size
%   maxfe       - max function evaluations (passed to the problem instance)
%   run         - number of independent runs to perform
%   save_foder  - folder path used to store a cached 'original_<Pro>_<Alg>.mat'
%   varargin    - additional args forwarded to the attacked algorithm constructor
%
% Output:
%   origin - array of baseline metric values (one entry per run). For
%            single-objective problems this is typically Min_value; for
%            multi-objective problems this is typically an HV value. For
%            problems that produce per-run populations, intermediate objective
%            values may be saved by parsave for later inspection.
%
% Behavior notes:
% - If a cached result file exists in the same folder as save_foder it is
%   loaded and returned immediately (avoids re-running expensive experiments).
% - This function preserves the original experiment flow and does not change
%   algorithm parameters; only documentation/comments were added.

function origin = Origin(Pro, Data, A_attacked, n, maxfe, run, save_foder, varargin)

% Build algorithm handle from the provided name and construct it with any
% additional parameters passed via varargin.
attackedAlgHandle = str2func(A_attacked);
Alg = attackedAlgHandle(varargin{:});

% Cache filename: original_<Pro>_<A_attacked>.mat placed next to the
% provided save folder. If present, return cached 'origin' immediately.
file = fullfile(fileparts(save_foder), sprintf("original_%s_%s.mat", Pro, A_attacked));

if exist(file, 'file') == 2
    loaded = load(file, 'origin');
    origin = loaded.origin;
    return;
end

% Main per-problem experiment loops. Each branch preserves the previous
% behavior: create the Attack* problem instance, call Alg.Solve(pro), and
% compute/store the baseline metric for that run.
if Pro == "TSP"

    % TSP: run 'run' independent runs and record Min_value of final result
    for i = 1:run
        pro = AttackTSP('N', n, 'maxFE', maxfe, 'parameter', {Data, 0, 0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end
    save(file, 'origin');

elseif Pro == "KP"
    % 0-1 Knapsack: Data contains values and weights concatenated; split them
    Data_P = Data(1:length(Data)/2);
    Data_W = Data(length(Data)/2+1:end);

    for i = 1:run
        % Keep original constructor call but note that n was previously
        % transposed (n' in original). Preserve original behavior.
        pro = AttackKP('N', n', 'maxFE', maxfe, 'parameter', {Data_P, Data_W, 0, 0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end

    save(file, 'origin');

elseif Pro == "MOTSP"

    % Multi-objective TSP: record hypervolume (HV) and also save intermediate
    % objective values per run for later analysis via parsave.
    for i = 1:run
        pro = AttackMOTSP('N', n, 'maxFE', maxfe, 'parameter', {Data, 0, 0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end}, [size(Data{1},1), size(Data{1},1)]);
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file, 'origin');

elseif Pro == "MOKP"
    % Multi-objective Knapsack: reshape Data into P and W and compute HV
    Data = reshape(Data, 4, []);
    Data_P = Data(1:2, :);
    Data_W = Data(3:4, :);

    for i = 1:run
        pro = AttackMOKP('N', n, 'maxFE', maxfe, 'parameter', {Data_P, Data_W, 0, 0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end}, sum(Data_P, 2)');
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file, 'origin');

elseif Pro == "PO"

    % Portfolio optimization: use problem's CalMetric('HV', result) and save
    for i = 1:run
        pro = AttackPO('N', n, 'maxFE', maxfe, 'parameter', {Data, 0, 0});
        Alg.Solve(pro);
        origin(i) = pro.CalMetric('HV', Alg.result{end});
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file, 'origin');

elseif Pro == "OPF"
    % Optimal power flow / case 118: record Min_value of final result
    for i = 1:run
        pro = AttackCase118('N', n, 'maxFE', maxfe, 'parameter', {Data, 0, 0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end
    save(file, 'origin');

elseif Pro == "CRT"

    % CRT (custom) branch: record HV and save per-run objectives
    for i = 1:run
        pro = AttackCRT('N', n, 'maxFE', maxfe, 'parameter', {[], 0, 0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end}, [1, 2]);
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file, 'origin');
end
end


function parsave(objs, fa_file, i, Pro, A_attacked)
% PARSAVE - helper to save per-run objective arrays produced during multi-
% objective experiments. Saves 'objs' into Val_<Pro>_<Alg>_<i>.mat in the
% given folder.
    file2 = sprintf('Val_%s_%s_%d.mat', Pro, A_attacked, i);
    file2 = fullfile(fa_file, file2);
    save(file2, "objs");
end