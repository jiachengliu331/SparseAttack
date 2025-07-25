classdef AttackMOKP_round < PROBLEM
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
        P;	% Profit of each item according to each knapsack
        W;  % Weight of each item according to each knapsack
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
            obj.encoding = ones(1,obj.D);
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            PopDec = round(PopDec);
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
            PopDec = round(PopDec);
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
    end
end