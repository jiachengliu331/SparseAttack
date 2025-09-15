classdef AttackPO < PROBLEM
    % <multi> <real> <large/none> <expensive/none> <sparse/none>
    % The portfolio optimization problem
    % dataNo --- 1 --- Number of dataset
    
    %------------------------------- Reference --------------------------------
    % Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
    % multi-objective optimization problems with sparse optimal solutions via
    % unsupervised neural networks, IEEE Transactions on Cybernetics, 2021,
    % 51(6): 3115-3128.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
    % for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    % The datasets are the minutely closing prices of EUR/CHF taken from MT4
    % No.   Name	Instruments   Length
    % 1     1000	   1000         50
    % 2     5000       5000         50
    
    properties
        %properties(Access = private)
        Data;
        Yield;
        Risk;
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Load data
            [obj.Data,obj.delta,obj.H] = obj.ParameterSet([],0,0);
            obj.Yield = log(obj.Data(:,2:end)) - log(obj.Data(:,1:end-1));
            obj.Risk  = cov(obj.Yield');
            % Parameter setting
            obj.M = 2;
            obj.D = size(obj.Yield,1);
            obj.lower    = zeros(1,obj.D) - 1;
            obj.upper    = zeros(1,obj.D) + 1;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            obj.Yield = log(obj.Data(:,2:end)) - log(obj.Data(:,1:end-1));
            obj.Risk  = cov(obj.Yield');
            PopDec = PopDec./repmat(max(sum(abs(PopDec),2),1),1,size(PopDec,2));
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                PopObj(i,1) = PopDec(i,:)*obj.Risk*PopDec(i,:)';
                PopObj(i,2) = 1 - sum(PopDec(i,:)*obj.Yield);
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopObj = Population.objs;
            Draw([PopObj(:,1),1-PopObj(:,2)],{'Risk','Return',[]});
        end
        %% Perturb solutions multiple times
        function PopX = Perturb(obj,PopDec,N)
            if nargin < 3 || N~=1
                N = obj.H;
                if N ==0 || N == 1
                    PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                else
                    temp = obj.Data;
                    for i = 1:N
                        Datamin = min(obj.Data,[],'all');
                        Datamax = max(obj.Data,[],'all');
                        Datanorm = (obj.Data - Datamin)/(Datamax - Datamin);

                        rand_perindex = randi(numel(Datanorm));

                        Delta = -Datanorm(rand_perindex) + rand();
                        Datanorm(rand_perindex) = Datanorm(rand_perindex) + Delta;

                        obj.Data = Datanorm * (Datamax - Datamin) + Datamin;

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
                        Datamin = min(obj.Data,[],'all');
                        Datamax = max(obj.Data,[],'all');
                        Datanorm = (obj.Data - Datamin)/(Datamax - Datamin);

                        rand_perindex = randi(numel(Datanorm));
                        Delta = -Datanorm(rand_perindex) + rand();

                        Datanorm(rand_perindex) = Datanorm(rand_perindex) + Delta;
                        obj.Data = Datanorm * (Datamax - Datamin) + Datamin;
                    end
                    PopX  = SOLUTION(obj.CalDec(PopDec),obj.CalObj(PopDec),obj.CalCon(PopDec));
                    obj.Data = temp;
                end
            end
        end    
    end
end