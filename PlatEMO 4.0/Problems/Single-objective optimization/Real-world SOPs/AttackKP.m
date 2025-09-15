classdef AttackKP < PROBLEM
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
        P;	        % Profit of each item
        W;          % Weight of each item
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [obj.P,obj.W,obj.delta,obj.H] = obj.ParameterSet([],[],0.01,100);
            obj.M = 1;
            obj.D = size(obj.P,2);
            obj.encoding = 4 + zeros(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = repmat(sum(obj.P,2)',size(PopDec,1),1) - PopDec*obj.P';
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = PopDec*obj.W' - sum(obj.W)/2;
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            [~,best] = min(Population.objs);
            Draw(obj.R(Population(best).dec([1:end,1]),:),'-k','LineWidth',1.5);
            Draw(obj.R);
        end
        %% Perturb solutions multiple times
        function PopX = Perturb(obj,PopDec)
            N = obj.H;
            if N == 0 || N == 1
                PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
            else
                temp1 = obj.P;
                temp2 = obj.W;
                for i = 1:N
                    Delta1 = rand(size(obj.P,1),size(obj.P,2));
                    Delta2 = rand(size(obj.W,1),size(obj.W,2));
                    obj.P = obj.P + Delta1;
                    obj.W = obj.W + Delta2;
                    PopX(i)  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                end
                obj.P = temp1;
                obj.W = temp2;
            end
        end
    end
end