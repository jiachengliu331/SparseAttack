function Population = FinalSelection(Population,Problem)
% Select final robust solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the convergence performance of solutions
    [PopObjV,PopConV]    = MeanEffective(Problem,Population);
    [~,rank]   = sort(FitnessSingleEEASRA(PopObjV,PopConV));
    Population = Population(rank(1:Problem.N));

end