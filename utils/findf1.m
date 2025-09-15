function PopDec_f1 = findf1(Population)
    PopObj = Population.objs;
    PopDec = Population.decs;
    N = size(PopObj,1);
    [FrontNo,~] = NDSort(PopObj,N);
    F1index = find(FrontNo==1);
    PopDec_f1 = PopDec(F1index,:);
end

