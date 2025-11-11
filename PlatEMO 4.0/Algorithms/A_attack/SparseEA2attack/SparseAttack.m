classdef SparseAttack < ALGORITHM
    % <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
    % SparseAttack
    
    %------------------------------- Reference --------------------------------
    % Y. Zhang, Y. Tian, and X. Zhang, Improved SparseEA for sparse large-scale
    % multi-objective optimization problems, Complex & Intelligent Systems,
    % 2021.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            [beta,n_s,t_s] = Algorithm.ParameterSet(0,0,0);

            %% Population initialization
            N = Problem.N;
            maxgen = ceil(Problem.maxFE/Problem.N);
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            for i = 1 : 5
               if i == 1
                    Dec        = repmat(Problem.upper,Problem.D,1);
                    Mask       = eye(Problem.D);
                    Population = Problem.Evaluation(Dec.*Mask);
                    TDec       = [TDec;Dec];
                    TMask      = [TMask;Mask];
                    TempPop    = [TempPop,Population];
                elseif i ==2
                    Dec        = repmat(Problem.lower,Problem.D,1);
                    Mask       = eye(Problem.D);
                    Population = Problem.Evaluation(Dec.*Mask);
                    TDec       = [TDec;Dec];
                    TMask      = [TMask;Mask];
                    TempPop    = [TempPop,Population];
                elseif i ==3
                    Dec        = repmat((Problem.upper - Problem.lower)/2,Problem.D,1);
                    Mask       = eye(Problem.D);
                    Population = Problem.Evaluation(Dec.*Mask);
                    TDec       = [TDec;Dec];
                    TMask      = [TMask;Mask];
                    TempPop    = [TempPop,Population];
                else
                    Dec        = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                    Dec(:,Problem.encoding==4) = 1;
                    Mask       = eye(Problem.D);
                    Population = Problem.Evaluation(Dec.*Mask);
                    TDec       = [TDec;Dec];
                    TMask      = [TMask;Mask];
                    TempPop    = [TempPop,Population];
                end
                Fitness        = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            % Generate initial population
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = false(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            indexno0_1 = find(sum(Mask==0,2)~=size(Mask,2));
            Dec = Dec(indexno0_1,:);
            Mask = Mask(indexno0_1,:);
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelectionAttack([Population],[Dec],[Mask],Problem.N,beta,t_s);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen              = ceil((Problem.FE)/N);
                K                = gen/maxgen;
                G                = max(floor((1-K)*N),floor(N/6));
                
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness);

                objs             = Population.objs;
                if nnz(objs(:,1) > t_s) < n_s
                    Offspring    = Problem.Evaluation(OffDec.*OffMask);
                else
                    F1index      = find(FrontNo==1);
                    F1mask       = Mask(F1index,:);

                    MatchOffIndex = find(ismember(OffMask,F1mask,'rows'));
                    MatchOffIndex = unique(MatchOffIndex);

                    Offfitness   = OffMask * Fitness';
                    index        = TournamentSelection(2,G,Offfitness');
                    index        = unique(index');
                    otherN       = N - length(index);

                    if length(MatchOffIndex) >= 0.5*otherN
                        index2 = TournamentSelection(2,floor(0.5*otherN),Offfitness(MatchOffIndex)');
                        MatchOffIndex = MatchOffIndex(index2');
                        MatchOffIndex = unique(MatchOffIndex);
                    end

                    indexoff = unique([MatchOffIndex;index]);
                    OffDec   = OffDec(indexoff,:);
                    OffMask  = OffMask(indexoff,:);

                    row_to_keep = any(OffMask~=0,2);
                    OffDec = OffDec(row_to_keep,:);
                    OffMask = OffMask(row_to_keep,:);

                    Offspring        = Problem.Evaluation(OffDec.*OffMask);
                end
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelectionAttack([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,beta,t_s);
            end
        end
    end
end