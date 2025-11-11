classdef AttackCase118 < PROBLEM
    % <many> <real> <large/none> <constrained/none>
    % Benchmark MOP with bias feature
    
    %------------------------------- Reference --------------------------------
    % Adaptive constraint differential evolution for optimal power flow
    % Shuijia Li , Wenyin Gong , Chengyu Hu a, Xuesong Yan a, Ling Wang b,
    % Qiong Gu
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    properties(Access = private)

        true_gencost;

        gencost;
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
        theta ;     % Sparsity of the Pareto set
        data;       %case118
    end
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.gencost,obj.delta,obj.H] = obj.ParameterSet([],0.01,100);
            %             obj.gencost = [obj.gencost,zeros(size(obj.gencost,1),1)];
            obj.M = 1;
            %             obj.M = 5;
            if isempty(obj.D); obj.D = 130; end
            %PG,Vg,Qc,T
            obj.lower    = [0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9];
            obj.upper    = [100	100	100	100	550	185	100	100	100	100	320	414	100	107	100	100	100	100	100	119	304	148	100	100	255	260	100	491	492	100	100	100	100	100	100	577	100	104	707	100	100	100	100	352	140	100	100	100	100	136	100	100	100 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 5 5 5 5 5 5 5 5 5 5 5 5 5 5 1.1 1.1 1.1 1.1 1.1 1.1 1.1 1.1 1.1];
            obj.encoding = ones(1,obj.D);
            obj.theta = cell(2,obj.N);
            obj.data = loadcase(case118);

            obj.true_gencost = obj.data.gencost(:,5:6);

            if ~isempty(obj.gencost)
                obj.data.gencost(:,5:6) = obj.gencost;
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Qbus = [5 34 37 44 45 46 48 74 79 82 83 105 107 110];%Number Q
            Tbranch = [8 32 36 51 93 95 102 107 127];            %Number T
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                obj.data.gen(1:29,2) = PopDec(i,1:29);
                obj.data.gen(31:54,2)= PopDec(i,31:54);
                obj.data.gen(1:54,6) = PopDec(i,54:107);
                obj.data.bus(Qbus,6) = PopDec(i,108:121);
                obj.data.branch(Tbranch,9) = PopDec(i,122:obj.D);
                mpopt = mpoption('verbose',0,'out.all',0);
                result = runpf(obj.data,mpopt);
                obj.theta{1,i} = result;
                obj.theta{2,i} = obj.data;
                %%TFC  total fuel cost objective
                PG = [PopDec(i,1:29),result.gen(30,2),PopDec(i,31:54)];
                costcoeff = obj.data.gencost(:,5:7);
                fuelcost = sum(costcoeff(:,1)+costcoeff(:,2).*PG'+costcoeff(:,3).*(PG.^2)');
                PopObj(i)= fuelcost/280000;
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            %             obj.data = loadcase(case118);
            %             Qbus = [5 34 37 44 45 46 48 74 79 82 83 105 107 110];%Number Q
            %             Tbranch = [8 32 36 51 93 95 102 107 127];              %Number T
            for i = 1 : size(PopDec,1)
                %             obj.data.gen(1:29,2) = PopDec(i,1:29);
                %             obj.data.gen(31:54,2)= PopDec(i,31:54);
                %             obj.data.gen(1:54,6) = PopDec(i,54:107);
                %             obj.data.bus(Qbus,6) = PopDec(i,108:121);
                %             obj.data.branch(Tbranch,9) = PopDec(i,122:obj.D);
                %             mpopt = mpoption('verbose',0,'out.all',0);
                %             result = runpf(obj.data,mpopt);
                result=obj.theta{1,i};
                obj.data=obj.theta{2,i};
                Vmax = obj.data.bus(:,12);
                Vmin = obj.data.bus(:,13);
                genbus = obj.data.gen(:,1);
                Qmax = obj.data.gen(:,4)/obj.data.baseMVA;
                Qmin = obj.data.gen(:,5)/obj.data.baseMVA;
                QG = result.gen(:,3)/obj.data.baseMVA;
                PGSmax = obj.data.gen(30,9);
                PGSmin = obj.data.gen(30,10);
                PGS = result.gen(30,2);
                %Generator 1 constraints
                PGS_err = (PGS<PGSmin)*(abs(PGSmin-PGS)/(PGSmax-PGSmin))+(PGS>PGSmax)*(abs(PGSmax-PGS)/(PGSmax-PGSmin));
                %Security constraints
                %             blimit = obj.data.branch(:,6);
                %             Slimit = sqrt(result.branch(:,14).^2+result.branch(:,15).^2);
                %             S_err = sum((Slimit>blimit).*abs(blimit-Slimit))/obj.data.baseMVA;
                %Q constraints
                %Q_err = (sum((QG<Qmin).*(abs(Qmin-QG)./(Qmax-Qmin))+(QG>Qmax).*(abs(Qmax-QG)./(Qmax-Qmin))))/118;
                %v constraints
                VI = result.bus(:,8);
                VI(genbus)=[];
                Vmax(genbus)=[];
                Vmin(genbus)=[];
                VI_err = (sum((VI<Vmin).*(abs(Vmin-VI)./(Vmax-Vmin))+(VI>Vmax).*(abs(Vmax-VI)./(Vmax-Vmin))))/100000;
                PopCon(i,:) = [VI_err,PGS_err];
            end
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [1,1,1,1,1];
        end
         %% Perturb solutions multiple times
        function PopX = Perturb(obj,PopDec,N)

            % PERTURB - generate perturbed evaluations for the OPF case (case118)
            %
            % This method creates N perturbed problem instances by modifying
            % the generator cost coefficients stored in obj.data.gencost(:,5:6),
            % evaluates the provided decision(s) `PopDec` on each perturbed
            % instance, and returns an array of SOLUTION objects. The original
            % implementation normalizes the true gencost, replaces a single
            % linear element with a uniform random value in [0,1), denormalizes
            % back to the original scale, and uses that modified gencost for
            % evaluation. That behavior is preserved here; only comments were
            % added to clarify the steps.

            % Use configured number of perturbation samples (H)
            N = obj.H;

            % If zero or one sample requested, evaluate with current data
            if N == 0 || N == 1
                PopX = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
            else
                % Save original gencost columns so we can restore them later
                temp = obj.data.gencost(:, 5:6);

                for i = 1:N
                    % Start from the stored true gencost (baseline)
                    tempgencost = obj.true_gencost;

                    % Compute min/max per column for normalization
                    max_1 = max(tempgencost(:, 1));
                    min_1 = min(tempgencost(:, 1));
                    max_2 = max(tempgencost(:, 2));
                    min_2 = min(tempgencost(:, 2));

                    % Normalize each column to [0,1]
                    tempgencostnorm = tempgencost;
                    tempgencostnorm(:, 1) = (tempgencost(:, 1) - min_1) / (max_1 - min_1);
                    tempgencostnorm(:, 2) = (tempgencost(:, 2) - min_2) / (max_2 - min_2);

                    % Choose a random linear index into the normalized matrix
                    % and replace that element with a random draw in [0,1).
                    % Note: this picks a single scalar entry (column-major
                    % linear index) rather than an entire row.
                    rand_perindex = randi(length(tempgencostnorm));
                    Delta = -tempgencostnorm(rand_perindex) + rand();
                    tempgencostnorm(rand_perindex) = tempgencostnorm(rand_perindex) + Delta;

                    % Denormalize back to the original scale
                    tempgencost(:, 1) = tempgencostnorm(:, 1) * (max_1 - min_1) + min_1;
                    tempgencost(:, 2) = tempgencostnorm(:, 2) * (max_2 - min_2) + min_2;

                    % Apply perturbed gencost to problem data and evaluate
                    obj.data.gencost(:, 5:6) = tempgencost;
                    PopX(i) = SOLUTION(obj.CalDec(PopDec), obj.CalObj(PopDec), obj.CalCon(PopDec));
                end

                % Restore original gencost columns
                obj.data.gencost(:, 5:6) = temp;
            end
        end
    end
end