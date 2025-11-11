classdef AttackMOTSP < PROBLEM
    % <multi/many> <permutation> <large/none>
    % The multi-objective traveling salesman problem
    % c --- 0 --- Correlation parameter
    
    %------------------------------- Reference --------------------------------
    % D. Corne and J. Knowles, Techniques for highly multiobjective
    % optimisation: some nondominated points are better than others,
    % Proceedings of the Annual Conference on Genetic and Evolutionary
    % Computation, 2007, 773-780.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [eOperator_spducational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    properties(Access = private)
        C;  % Adjacency matrix of each map
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
        Data;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [obj.C,obj.delta,obj.H] = obj.ParameterSet([],0.05,50);
            obj.M = length(obj.C);
            obj.D = size(obj.C{1},1);
            obj.encoding = 5 + zeros(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            PopObj = zeros(N,obj.M);
            for i = 1 : obj.M
                for j = 1 : N
                    for k = 1 : D-1
                        PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,k),PopDec(j,k+1));
                    end
                    PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,D),PopDec(j,1));
                end
            end
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = zeros(1,obj.M) + obj.D;
        end
        %% Perturb solutions multiple times
        %
        % This method perturbs the distance/cost matrices (stored in obj.C)
        % and computes perturbed solutions. It intentionally supports two
        % related modes (kept from the original implementation):
        %
        % Mode A (sample-per-perturb): when called without an explicit N or
        % when the input condition (nargin < 3 || N ~= 1) is true, the code
        % sets N = obj.H and for each sample it temporarily perturbs one
        % matrix element, evaluates the solution, and then restores the
        % original data. This produces N independently perturbed solutions
        % stored in PopX (one per perturbation).
        %
        % Mode B (cumulative-perturb then evaluate once): when the alternate
        % branch is taken the code applies N random perturbations cumulatively
        % to the matrices, evaluates once after all changes, then restores
        % the original matrices. This yields a single evaluation that
        % reflects the compound effect of multiple small changes.
        %
        % Implementation notes:
        % - obj.Data is set to obj.C (the per-map adjacency matrices). The
        %   matrices are stored in a compact upper-triangular vector form
        %   for the chosen indexing; the code uses tril/triu mirror to keep
        %   the matrices symmetric after modifying the compact storage.
        % - len_d is the number of entries in the stored triangular vector
        %   representation (used to pick a random index to perturb).
        % - Delta is computed as (-x + rand()) so the updated element
        %   becomes rand() (i.e., a uniform random value in (0,1)). This
        %   matches the original behavior where a single element is replaced
        %   by rand().
        % - temp holds the original Data so we can restore it after
        %   temporary perturbations and avoid permanent mutation of the
        %   problem instance.
        function PopX = Perturb(obj,PopDec,N)
            % Begin by copying adjacency matrices into a working field
            obj.Data = obj.C;
            % Number of compacted (lower-triangular) elements in one map
            len_d = length(find(tril(obj.Data{1}, -1) ~= 0));

            % Preserve the original branching logic; the outer condition
            % controls which of the two modes is used.
            if nargin < 3 || N ~= 1
                % Mode A entry: sample-per-perturb
                N = obj.H;
                if N == 0 || N == 1
                    % No perturbation requested -> evaluate once
                    PopX = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                else
                    % Save original data for restoration after each sample
                    temp = obj.Data;
                    for i = 1:N
                        rand_m = randi(obj.M);    % which map to perturb
                        rand_d = randi(len_d);    % which compact-index to replace

                        % Compute Delta so that obj.Data{rand_m}(rand_d) becomes rand()
                        Delta = -obj.Data{rand_m}(rand_d) + rand();

                        % Apply the replacement and rebuild the symmetric matrix
                        obj.Data{rand_m}(rand_d) = obj.Data{rand_m}(rand_d) + Delta;
                        obj.Data{rand_m} = tril(obj.Data{rand_m}, -1) + triu(obj.Data{rand_m}', 1);

                        % Evaluate the solution under this temporary change
                        PopX(i) = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                        % Restore original data before next sample
                        obj.Data = temp;
                    end
                end
            else
                % Mode B entry: cumulative perturbations then single evaluate
                if N == 0 || N == 1
                    PopX = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                else
                    temp = obj.Data;
                    for i = 1:N
                        rand_m = randi(obj.M);
                        rand_d = randi(len_d);

                        Delta = -obj.Data{rand_m}(rand_d) + rand();

                        obj.Data{rand_m}(rand_d) = obj.Data{rand_m}(rand_d) + Delta;
                        obj.Data{rand_m} = tril(obj.Data{rand_m}, -1) + triu(obj.Data{rand_m}', 1);
                    end
                    % One evaluation after cumulative changes
                    PopX = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                    % Restore the original matrices
                    obj.Data = temp;
                end
            end
        end
    end
end