classdef AttackMOKP < PROBLEM
    % <multi/many> <binary> <large/none>
    % The multi-objective knapsack problem
    
    %------------------------------- Reference --------------------------------
    % E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
    % comparative case study and the strength Pareto approach, IEEE
    % Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    properties(Access = private)
        P;	        % Profit of each item according to each knapsack
        W;          % Weight of each item according to each knapsack
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [obj.P,obj.W,obj.delta,obj.H] = obj.ParameterSet([],[],0.05,50);
            obj.M = size(obj.P,1);
            obj.D = size(obj.P,2);
            obj.encoding = 3 + ones(1,obj.D);
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            C = sum(obj.W,2)/2;
            [~,rank] = sort(max(obj.P./obj.W));
            for i = 1 : size(PopDec,1)
                while any(obj.W*PopDec(i,:)'>C)
                    k = find(PopDec(i,rank),1);
                    PopDec(i,rank(k)) = 0;
                end
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = repmat(sum(obj.P,2)',size(PopDec,1),1) - PopDec*obj.P';
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = sum(obj.P,2)';
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.decs*obj.P',{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
        end
        %% Perturb solutions multiple times
        %
        % This method supports two related modes (preserved from the original
        % implementation):
        %
        % Mode A (sample-per-perturb): when called without an explicit N or
        % when the input test (nargin < 3 || N ~= 1) is true, the code sets
        % N = obj.H and produces N perturbed solutions. For each sample the
        % code temporarily perturbs either a single profit (P) or weight (W)
        % entry, computes the solution (via CalDec/CalObj/CalCon), then
        % restores the original P/W before generating the next sample. This
        % yields an array PopX of length N where each element was generated
        % by exactly one small, temporary perturbation.
        %
        % Mode B (cumulative-perturb then evaluate once): when the initial
        % condition is false (the 'else' branch in the original code), the
        % method again sets N = obj.H but instead applies N random
        % perturbations cumulatively to P and/or W, computes one solution
        % from the cumulatively perturbed problem, and then restores P and
        % W. This yields a single-element PopX (one evaluation after many
        % small changes) and is different from Mode A.
        %
        % Notes:
        %  - The two branches intentionally differ in semantics. I preserved
        %    that behavior exactly but added comments to explain what's
        %    happening. If you prefer a single clear policy (recommended),
        %    we can refactor to a single-mode implementation.
        %  - temp1/temp2 are used to save and restore the original P and W
        %    values so that perturbations do not permanently mutate the
        %    problem instance.
        %  - Delta is computed as: Delta = -x + rand() which produces a
        %    random value in ( -x, 1-x ). The code then adds Delta to the
        %    selected P/W element, effectively replacing it with a random
        %    number in (0,1) (since new value = x + (-x + rand()) = rand()).
        %    In short: the perturbation sets the chosen entry to rand().
        function PopX = Perturb(obj,PopDec,N)
            % If N is not provided or N~=1 the original code sets sampling
            % behavior handled in the first branch. Preserve that check.
            if nargin < 3 || N~=1
                N = obj.H;
                if N == 0 || N == 1
                    % No perturbation requested: just evaluate once and wrap
                    % the result in a SOLUTION object.
                    PopX  = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                else
                    % Mode A: sample-per-perturb. Save originals so we can
                    % restore after each temporary perturbation.
                    temp1 = obj.P;
                    temp2 = obj.W;
                    for i = 1:N
                        rand_m = randi(2);           % choose whether to perturb P or W
                        if rand_m == 1
                            rand_d = randi(obj.D);   % which dimension/item
                            % compute Delta so that obj.P(rand_d) becomes rand()
                            Delta = -obj.P(rand_d) + rand();
                            obj.P(rand_d) = obj.P(rand_d) + Delta;
                            % evaluate with this temporary change
                            PopX(i)  = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                            % restore original profits
                            obj.P = temp1;
                        else
                            rand_d = randi(obj.D);
                            Delta = -obj.W(rand_d) + rand();
                            obj.W(rand_d) = obj.W(rand_d) + Delta;
                            PopX(i)  = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                            % restore original weights
                            obj.W = temp2;
                        end
                    end
                end
            else
                % The else branch preserves the original alternate behavior
                % (Mode B): apply N perturbations cumulatively, evaluate
                % once, then restore the original P/W.
                N = obj.H;
                if N == 0 || N == 1
                    PopX  = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                else
                    % Save originals
                    temp1 = obj.P;
                    temp2 = obj.W;
                    % Apply N small random replacements to entries of P or W
                    for i = 1:N
                        rand_m = randi(2);
                        if rand_m == 1
                            rand_d = randi(obj.D);
                            Delta = -obj.P(rand_d) + rand();
                            obj.P(rand_d) = obj.P(rand_d) + Delta;
                        else
                            rand_d = randi(obj.D);
                            Delta = -obj.W(rand_d) + rand();
                            obj.W(rand_d) = obj.W(rand_d) + Delta;
                        end
                    end
                    % One evaluation after cumulative changes
                    PopX  = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                    % Restore originals
                    obj.P = temp1;
                    obj.W = temp2;
                end
            end
        end
    end
end