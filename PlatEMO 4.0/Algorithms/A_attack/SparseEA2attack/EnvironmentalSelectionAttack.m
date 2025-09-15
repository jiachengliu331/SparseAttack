function [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelectionAttack(Population,Dec,Mask,N,beta,t_s)
% The environmental selection of SparseEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicated solutions
    [~,uni]    = unique(Population.objs,'rows');
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    objs       = Population.objs;

    indexno    = find(objs(:,1)~=0);
    Population = Population(indexno);
    objs       = Population.objs;
    Dec        = Dec(indexno,:);
    Mask       = Mask(indexno,:);

    objs(:,2)  = objs(:,2) + beta*size(Dec,2)*(max(t_s,objs(:,1)) - t_s).*abs(objs(:,2));
    N          = min(N,length(Population));

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(objs,Population.cons,N);
    Next       = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis   = CrowdingDistance(objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last       = find(FrontNo==MaxFNo);
    [~,Rank]   = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
end