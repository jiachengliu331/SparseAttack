classdef SparseAttack < ALGORITHM
    % <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
    % SparseAttack

    %--------------------------------------------------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference:
    % Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum],
    % IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87.
    %--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            [beta, n_s, t_s] = Algorithm.ParameterSet(0, 0, 0);

            %% Initialization
            N = Problem.N;
            maxgen = ceil(Problem.maxFE / N);

            TDec    = [];   % Temporary decision variables
            TMask   = [];   % Temporary masks
            TempPop = [];   % Temporary population
            Fitness = zeros(1, Problem.D);  % Fitness per dimension

            % Generate reference individuals for each extreme decision setting
            for i = 1 : 5
                if i == 1
                    % All variables at upper bound
                    Dec = repmat(Problem.upper, Problem.D, 1);
                elseif i == 2
                    % All variables at lower bound
                    Dec = repmat(Problem.lower, Problem.D, 1);
                elseif i == 3
                    % All variables at midpoint
                    Dec = repmat((Problem.upper - Problem.lower) / 2, Problem.D, 1);
                else
                    % Random values within bounds
                    Dec = unifrnd(repmat(Problem.lower, Problem.D, 1), ...
                                  repmat(Problem.upper, Problem.D, 1));
                end

                Mask = eye(Problem.D);  % One variable active per solution
                Population = Problem.Evaluation(Dec .* Mask);  % Evaluation

                TDec    = [TDec; Dec];
                TMask   = [TMask; Mask];
                TempPop = [TempPop, Population];

                % Update dimension-wise fitness using non-dominated sorting
                Fitness = Fitness + NDSort([Population.objs, Population.cons], inf);
            end

            %% Initial Population
            % Generate random decision values within bounds
            Dec = unifrnd(repmat(Problem.lower, N, 1), repmat(Problem.upper, N, 1));
            Dec(:, Problem.encoding == 4) = 1;

            % Randomly generate sparse masks with selected active variables
            Mask = false(N, Problem.D);
            for i = 1 : N
                % Random number of active variables (via tournament selection)
	Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end

            % Remove individuals with all-zero masks (no active variables)
            indexValid = find(sum(Mask == 0, 2) ~= size(Mask, 2));
            Dec  = Dec(indexValid, :);
            Mask = Mask(indexValid, :);

            % Evaluate initial sparse population
            Population = Problem.Evaluation(Dec .* Mask);

            % First environmental selection
            [Population, Dec, Mask, FrontNo, CrowdDis] = ...
                EnvironmentalSelectionAttack([Population, TempPop], ...
                                             [Dec; TDec], [Mask; TMask], ...
                                             N, beta, t_s);

            %% Optimization loop
            while Algorithm.NotTerminated(Population)
                gen = ceil(Problem.FE / N);     % Current generation
                K   = gen / maxgen;             % Progress ratio
                G   = max(floor((1 - K) * N), floor(N / 6)); 

                % Parent selection
                MatingPool = TournamentSelection(2, 2 * N, FrontNo, -CrowdDis);
                [OffDec, OffMask] = Operator(Problem, Dec(MatingPool, :), ...
                                             Mask(MatingPool, :), Fitness);

                % If too few solutions exceed threshold, evaluate all offspring
                objs = Population.objs;
                if nnz(objs(:, 1) > t_s) < n_s
                    Offspring = Problem.Evaluation(OffDec .* OffMask);
                else
                    F1index = find(FrontNo == 1);    % First front
                    F1mask  = Mask(F1index, :);      % Masks of top individuals

                    MatchOffIndex = find(ismember(F1mask, OffMask, 'rows'));
                    MatchOffIndex = unique(MatchOffIndex);

                    Offfitness = OffMask * Fitness';

                    index = TournamentSelection(2, G, Offfitness');
                    index = unique(index');
                    otherN = N - length(index);

                    if length(MatchOffIndex) >= 0.5 * otherN
                        index2 = TournamentSelection(2, floor(0.5 * otherN), Offfitness);
                        MatchOffIndex = unique(index2');
                    end

                    indexoff = unique([MatchOffIndex; index]);
                    OffDec   = OffDec(indexoff, :);
                    OffMask  = OffMask(indexoff, :);

                    % Remove masks with no active variables
                    row_to_keep = any(OffMask ~= 0, 2);
                    OffDec  = OffDec(row_to_keep, :);
                    OffMask = OffMask(row_to_keep, :);

                    % Evaluate selected offspring
                    Offspring = Problem.Evaluation(OffDec .* OffMask);
                end

                % Environmental selection for next generation
                [Population, Dec, Mask, FrontNo, CrowdDis] = ...
                    EnvironmentalSelectionAttack([Population, Offspring], ...
                                                 [Dec; OffDec], ...
                                                 [Mask; OffMask], ...
                                                 N, beta, t_s);
            end
        end
    end
end
