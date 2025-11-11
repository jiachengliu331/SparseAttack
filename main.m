function main(varargin)

    % MAIN - entry point for running sparse-attack experiments using PlatEMO
    % This script parses optional parameters, loads dataset/problem data,
    % constructs objective functions (one measuring sparsity and another
    % the perturbation objective), runs the attack algorithm repeatedly
    % via PlatEMO, and saves the resulting decision variables and
    % objective values to disk.

    % ------------------------
    % Parse optional inputs
    % ------------------------
    parser = inputParser;

    % General experiment identifiers
    addOptional(parser, "maxnumber", 1);

    % Problem and algorithm selection (default values chosen here)
    % Note: using @ischar restricts inputs to character arrays; if you
    % prefer MATLAB string objects ("..."), consider allowing isstring.
    addOptional(parser, "Pro", 'TSP', @ischar);
    addOptional(parser, "A_attacked", 'GA', @ischar);
    addOptional(parser, "A_attack", 'SparseAttack', @ischar);

    % Resource limits and population sizes
    addOptional(parser, "maxfe", 2500, @(x) x > 0);   % max function evaluations for attacked algorithm
    addOptional(parser, "MAXFE1", 3000, @(x) x > 0);  % max function evaluations for attack algorithm
    % Population sizes
    addOptional(parser, "n", 30, @(x) x > 0);   % population size of the attacked algorithm
    addOptional(parser, "N1", 30, @(x) x > 0);  % population size of the attack algorithm

    % SparseAttack hyperparameters 
    addOptional(parser, "beta", 0.05);
    addOptional(parser, "n_s", 3);
    addOptional(parser, "t_s", 0.1);

    % Number of independent runs
    addOptional(parser, "run", 10);

    % Perturbation scale
    addOptional(parser, "scalefactor", 1);

    parse(parser, varargin{:});

    % ------------------------
    % Build config struct
    % ------------------------
    config.maxnumber   = parser.Results.maxnumber;

    config.Pro         = parser.Results.Pro;
    config.A_attacked  = parser.Results.A_attacked;
    config.A_attack    = parser.Results.A_attack;

    config.maxfe       = parser.Results.maxfe;
    config.MAXFE       = parser.Results.MAXFE1;
    config.n           = parser.Results.n;
    config.N           = parser.Results.N1;

    config.beta        = parser.Results.beta;
    config.n_s         = parser.Results.n_s;

    config.t_s         = parser.Results.t_s;

    config.run         = parser.Results.run;

    config.scalefactor = parser.Results.scalefactor;

    % ------------------------
    % Setup paths and save folder
    % ------------------------
    mainfolder = fileparts(mfilename("fullpath"));
    % Add entire project tree to MATLAB path so helper functions are found
    addpath(genpath(mainfolder));

    % Folder structure: <mainfolder>/<Pro>/test<maxnumber>/<A_attacked>_<A_attack>/
    save_folder_pro = fullfile(mainfolder, config.Pro);
    save_folder     = fullfile(save_folder_pro, sprintf("test%d", config.maxnumber));
    save_folder     = fullfile(save_folder, sprintf("%s_%s", config.A_attacked, config.A_attack));
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Persist the configuration used for this experiment
    save(fullfile(save_folder, "config.mat"), '-struct', 'config');

    % ------------------------
    % Load data and prepare objectives
    % ------------------------
    % loadOrgData should return problem Data and encoding/lower/upper bounds
    [Data, Encoding, Lower, Upper] = loadOrgData(config);

    % Origin computes a baseline (e.g., original solution/objectives) for
    % the attacked algorithm. The call returns a vector; we take the mean
    % as done originally.
    origin = Origin(config.Pro, Data, config.A_attacked, config.n, config.maxfe, config.run, save_folder, 'save', 0);
    origin = mean(origin);

    % Fun1: sparsity objective — mean proportion of non-zero decision vars
    Fun1 = @(x) mean(x ~= 0);

    % Fun2: perturbation objective constructed by Perturbation().
    Fun2 = Perturbation(config.Pro, Data, config.A_attacked, origin, config.n, config.maxfe);

    % ------------------------
    % Run the attack algorithm multiple times and save results
    % ------------------------
    parfor fileno = 1:config.run
        % Call platemo with the attack algorithm and its parameters.
        % 'algorithm' takes a cell where the first element is the function
        % handle to the algorithm and the remaining elements are algorithm-specific params.
        [DEC_attack, OBJ_attack, ~, ~] = platemo('algorithm', {str2func(config.A_attack), config.beta, config.n_s, config.t_s}, ...
            'N', config.N, 'maxFE', config.MAXFE, 'save', 10, ...
            'objFcn', {Fun1, Fun2}, 'encoding', Encoding, 'lower', Lower, 'upper', Upper);

        % Persist the decision variables and objective values for this run
        parsave(DEC_attack, OBJ_attack, config.A_attacked, config.A_attack, fileno, save_folder);
    end

end

function parsave(DEC_attack, OBJ_attack, AttackedAlgorithm, AttackAlgorithm, fileno, save_folder)
    % PERSAVE - helper to save DEC and OBJ results for a single run
    % Files written:
    %   DEC_attack_<AttackedAlgorithm>_<AttackAlgorithm>_<fileno>.mat
    %   OBJ_attack_<AttackedAlgorithm>_<AttackAlgorithm>_<fileno>.mat

    file2 = sprintf('DEC_attack_%s_%s_%d.mat', AttackedAlgorithm, AttackAlgorithm, fileno);
    file2 = fullfile(save_folder, file2);
    save(file2, 'DEC_attack');

    file3 = sprintf('OBJ_attack_%s_%s_%d.mat', AttackedAlgorithm, AttackAlgorithm, fileno);
    file3 = fullfile(save_folder, file3);
    save(file3, 'OBJ_attack');
end