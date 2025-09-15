function Func2 = Perturbation(Pro, Data, A_attacked, origin, N, MAXFE)
    
    if Pro == "TSP"
        
        Func2 = @(x)TSP_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "KP"
    
        Func2 = @(x)KP_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "MOTSP"
    
        Func2 = @(x)MOTSP_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "MOKP"
    
        Func2 = @(x)MOKP_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "PO"
    
        Func2 = @(x)PO_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "OPF"
    
        Func2 = @(x)OPF_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    elseif Pro == "CRT"
    
        Func2 = @(x)CRT_Perturbation(x, Data, A_attacked, origin, N, MAXFE);
    
    end

end


function [Diff,Obj2] = TSP_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);
    
    H = sum(x~=0);
    delta = 0;

    Pro = AttackTSP('N',N,'maxFE',MAXFE,'parameter',{Data+reshape(x,[],2),delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);

    Obj2 = calOriginObj_Single(Alg, Data, @TSP_calperobj);
    Diff = (-Obj2 + origin)/origin;

end

function [Diff,Obj2] = KP_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);

    H = sum(x~=0);
    delta = 0;

    Datatemp = Data+x;
    DataP_P = Datatemp(1:length(x)/2);
    DataP_W = Datatemp(length(x)/2+1:end);

    Pro = AttackKP('N',N','maxFE',MAXFE,'parameter',{DataP_P,DataP_W,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);

    Data_P = Data(1:length(x)/2);
    Obj2 = calOriginObj_Single(Alg, Data_P, @KP_calperobj);

    Diff = (-Obj2 + origin)/origin;

end

function [Diff,pop] = MOTSP_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);
    
    H = sum(x~=0);
    delta =0;

    index1 = find(tril(Data{1},-1)~=0);
    index2 = find(tril(Data{2},-1)~=0);
    
    Datap = Data;
    i=1;
    for j = 1:length(index1)
        Datap{1}(index1(j)) = Data{1}(index1(j)) + x(i);
        i  = i+1;
    end
    Datap{1} = tril(Datap{1},-1) + triu(Datap{1}',1);
    for k = 1:length(index2)
        Datap{2}(index2(k)) = Data{2}(index2(k)) + x(i);
        i = i+1;
    end
    Datap{2} = tril(Datap{2},-1) + triu(Datap{2}',1);
    
    Pro = AttackMOTSP('N',N,'maxFE',MAXFE,'parameter',{Datap,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);
    [HV2,pop] = calOriginObj_Multi(Alg, Data, @MOTSP_calperobj,[size(Data{1},1),size(Data{1},1)]);
    Diff = (HV2 - origin)/origin;

end

function [Diff,pop] = MOKP_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);

    H = sum(x~=0);
    delta =0;

    Data    = reshape(Data,4,[]);
    x       = reshape(x,4,[]);
    Data_P  = Data(1:2,:);
    Datap   = Data+x;
    DataP_P = Datap(1:2,:);
    DataP_W = Datap(3:4,:);

    Pro = AttackMOKP('N',N,'maxFE',MAXFE,'parameter',{DataP_P,DataP_W,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);
    
    [HV2,pop] = calOriginObj_Multi(Alg, Data, @MOKP_calperobj,sum(Data_P,2)');
    Diff = (HV2 - origin)/origin;

end

function [Diff,pop] = PO_Perturbation(x, Data, A_attacked, origin, N, MAXFE)
    
    attackedAlgHandle = str2func(A_attacked);
    
    H = sum(x~=0);
    delta =0;
    
    Datamin = min(Data,[],'all');
    Datamax = max(Data,[],'all');
    Datanorm = (Data - Datamin)/(Datamax - Datamin);
    Datanormp = Datanorm+reshape(x,[],size(Datanorm,2));
    Datap = Datanormp*(Datamax - Datamin) + Datamin;
    
    Pro = AttackPO('N',N,'maxFE',MAXFE,'parameter',{Datap,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);

    [HV2,pop] = calOriginObj_Multi(Alg, Data, @PO_calperobj,Pro.optimum);
    Diff = (HV2 - origin)/origin;

end

function [Diff,Obj2] = OPF_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);
    
    H = sum(x~=0);
    delta = 0;
    
    max_1      = max(Data(:,1));
    min_1      = min(Data(:,1));
    max_2      = max(Data(:,2));
    min_2      = min(Data(:,2));
    Datanorm(:,1) = (Data(:,1) - min_1)/(max_1 - min_1);
    Datanorm(:,2) = (Data(:,2) - min_2)/(max_2 - min_2);

    Datanormp = Datanorm + reshape(x,[],2);
    Datap(:,1) = Datanormp(:,1) *(max_1 - min_1) + min_1;
    Datap(:,2) = Datanormp(:,2) *(max_2 - min_2) + min_2;
    
    Pro = AttackCase118('N',N,'maxFE',MAXFE,'parameter',{Datap,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);

    Obj2 = calOriginObj_Single(Alg, Data, @OPF_calperobj);
    Diff = (-Obj2 + origin)/origin;
    
end

function [Diff,pop] = CRT_Perturbation(x, Data, A_attacked, origin, N, MAXFE)

    attackedAlgHandle = str2func(A_attacked);

    H = sum(x~=0);
    delta = 0;

    Datanorm   = (Data - min(Data))/(max(Data)-min(Data));
    Datanormp = Datanorm+x';
    Datap = Datanormp*(max(Data)-min(Data))+min(Data);

    Pro = AttackCRT('N',N,'maxFE',MAXFE,'parameter',{Datap,delta,H});
    Alg = attackedAlgHandle('save',0);
    Alg.Solve(Pro);

    [HV2,pop] = calOriginObj_Multi(Alg, Data, @CRT_calperobj,[1,2]);
    Diff = (HV2 - origin)/origin;

end


function PopObj = TSP_calperobj(PopDec,Data)

    C = pdist2(Data,Data);
    PopObj = zeros(size(PopDec,1),1);
    for i = 1 : size(PopDec,1)
        for j = 1 : size(PopDec,2)-1
            PopObj(i) = PopObj(i) + C(PopDec(i,j),PopDec(i,j+1));
        end
        PopObj(i) = PopObj(i) + C(PopDec(i,end),PopDec(i,1));
    end

end


function PopObj = KP_calperobj(PopDec,Data_P)

    P = Data_P;
    PopObj = repmat(sum(P,2)',size(PopDec,1),1) - PopDec*P';

end


function PopObj = MOTSP_calperobj(PopDec,Data)

    [N,D]  = size(PopDec);
    PopObj = zeros(N,length(Data));
    for i = 1 : length(Data)
        for j = 1 : N
            for k = 1 : D-1
                PopObj(j,i) = PopObj(j,i) + Data{i}(PopDec(j,k),PopDec(j,k+1));
            end
            PopObj(j,i) = PopObj(j,i) + Data{i}(PopDec(j,D),PopDec(j,1));
        end
    end

end


function PopObj = MOKP_calperobj(PopDec,Data)

    Data_P = Data(1:2,:);
    Data_W = Data(3:4,:);
    C = sum(Data_W,2)/2;
    [~,rank] = sort(max(Data_P./Data_W));
    for i = 1 : size(PopDec,1)
        while any(Data_W*PopDec(i,:)'>C)
            k = find(PopDec(i,rank),1);
            PopDec(i,rank(k)) = 0;
        end
    end
    PopObj = repmat(sum(Data_P,2)',size(PopDec,1),1) - PopDec*Data_P';

end


function PopObj = PO_calperobj(PopDec,Data)

    Yield = log(Data(:,2:end)) - log(Data(:,1:end-1));
    Risk  = cov(Yield');
    PopDec = PopDec./repmat(max(sum(abs(PopDec),2),1),1,size(PopDec,2));
    PopObj = zeros(size(PopDec,1),2);
    for i = 1 : size(PopDec,1)
        PopObj(i,1) = PopDec(i,:)*Risk*PopDec(i,:)';
        PopObj(i,2) = 1 - sum(PopDec(i,:)*Yield);
    end

end


function PopObj = OPF_calperobj(PopDec, Data)

    data = loadcase(case118);
    
    Qbus = [5 34 37 44 45 46 48 74 79 82 83 105 107 110];%Number Q
    Tbranch = [8 32 36 51 93 95 102 107 127];            %Number T
    PopObj = zeros(size(PopDec,1),1);
    for i = 1 : size(PopDec,1)
        data.gen(1:29,2) = PopDec(i,1:29);
        data.gen(31:54,2)= PopDec(i,31:54);
        data.gen(1:54,6) = PopDec(i,54:107);
        data.bus(Qbus,6) = PopDec(i,108:121);
        data.branch(Tbranch,9) = PopDec(i,122:130);
        mpopt = mpoption('verbose',0,'out.all',0);
        result = runpf(data,mpopt);
    
        PG = [PopDec(i,1:29),result.gen(30,2),PopDec(i,31:54)];
        costcoeff = data.gencost(:,5:7);
        fuelcost = sum(costcoeff(:,1)+costcoeff(:,2).*PG'+costcoeff(:,3).*(PG.^2)');
        PopObj(i)=fuelcost/280000;
    end
end



