classdef AttackTSP < PROBLEM
    % <single> <permutation> <large/none>
    % The traveling salesman problem
    
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
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    properties
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
        R;          % Locations of points
        C;          % Adjacency matrix
        Data;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [obj.R,obj.delta,obj.H] = obj.ParameterSet([],0.01,100);
            obj.M = 1;
            obj.D = size(obj.R,1);
            obj.encoding = 5 + zeros(1,obj.D);
            obj.C = pdist2(obj.R,obj.R);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            obj.C = pdist2(obj.R,obj.R);
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                for j = 1 : size(PopDec,2)-1
                    PopObj(i) = PopObj(i) + obj.C(PopDec(i,j),PopDec(i,j+1));
                end
                PopObj(i) = PopObj(i) + obj.C(PopDec(i,end),PopDec(i,1));
            end
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            [~,best] = min(Population.objs);
            Draw(obj.R(Population(best).dec([1:end,1]),:),'-k','LineWidth',1.5);
            Draw(obj.R);
        end
        function PopX = Perturb(obj, PopDec)

            % PERTURB - generate perturbed evaluations for given decision(s)
            %
            % This method produces one or more perturbed SOLUTION objects by
            % modifying the coordinate matrix `obj.R` and evaluating the
            % provided `PopDec`. The original implementation chooses a
            % random linear index in `obj.R` and replaces that element with
            % a uniform random value in [0,1). That behavior is preserved
            % here; only comments were added to document the logic.

            % If H is 0 or 1, return a single evaluation using the current R
            if obj.H == 0 || obj.H == 1
                PopX = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
            else
                % Preserve original city coordinates and restore after each
                % temporary perturbation to avoid permanently changing obj.R
                temp = obj.R;

                % Generate H perturbed variants
                for i = 1:obj.H
                    perturb_R = obj.R; % local copy

                    % Choose a random linear index into the coordinate matrix
                    % (note: this picks either an x or y coordinate of a city)
                    rand_perindex = randi(numel(obj.R));
                    Delta = -perturb_R(rand_perindex) + rand();
                    perturb_R(rand_perindex) = perturb_R(rand_perindex) + Delta;

                    % Temporarily use the perturbed coordinates and evaluate
                    obj.R = perturb_R;
                    PopX(i) = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));

                    % Restore original coordinates
                    obj.R = temp;
                end
            end

        end
    end
end