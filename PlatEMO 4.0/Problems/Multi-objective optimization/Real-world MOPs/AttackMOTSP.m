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
        function PopX = Perturb(obj,PopDec,N)
            obj.Data = obj.C;
            len_d = length(find(tril(obj.Data{1},-1)~=0));
            if nargin < 3 || N~=1
                N = obj.H;
                if N == 0 || N == 1
                    PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                else
                    temp = obj.Data;
                    for i = 1:N
                        rand_m = randi(obj.M);
                        rand_d = randi(len_d);

                        Delta = -obj.Data{rand_m}(rand_d) + rand();

                        obj.Data{rand_m}(rand_d) = obj.Data{rand_m}(rand_d) + Delta;
                        obj.Data{rand_m} = tril(obj.Data{rand_m},-1) + triu(obj.Data{rand_m}',1);
                        
                        PopX(i)  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                        obj.Data = temp;
                    end
                end
            else
                if N ==0 || N == 1
                    PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                else
                    temp = obj.Data;
                    for i = 1:N
                        rand_m = randi(obj.M);
                        rand_d = randi(len_d);

                        Delta = -obj.Data{rand_m}(rand_d) + rand();

                        obj.Data{rand_m}(rand_d) = obj.Data{rand_m}(rand_d) + Delta;
                        obj.Data{rand_m} = tril(obj.Data{rand_m},-1) + triu(obj.Data{rand_m}',1);
                    end
                    PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                    obj.Data = temp;
                end
            end
        end
    end
end