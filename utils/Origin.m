function origin = Origin(Pro, Data, A_attacked, n, maxfe, run, save_foder, varargin)

attackedAlgHandle = str2func(A_attacked);
Alg = attackedAlgHandle(varargin{:});

file = fullfile(fileparts(save_foder),sprintf("original_%s_%s.mat",Pro, A_attacked));

if exist(file,'file') == 2
    origin = load(file,'origin');
    origin = origin.origin;
    return;
end

if Pro == "TSP"

    for i=1:run
        pro = AttackTSP('N', n, 'maxFE', maxfe, 'parameter', {Data, 0, 0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end
    save(file,'origin');

elseif Pro == "KP"
    
    Data_P = Data(1:length(Data)/2);
    Data_W = Data(length(Data)/2+1:end);
    
    for i=1:run
        pro = AttackKP('N',n','maxFE',maxfe,'parameter',{Data_P,Data_W,0,0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end
    
    save(file,'origin');

elseif Pro == "MOTSP"

    for i=1:run
        pro = AttackMOTSP('N',n,'maxFE',maxfe,'parameter',{Data,0,0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end},[size(Data{1},1),size(Data{1},1)]);
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file,'origin');

elseif Pro == "MOKP"
    
    Data = reshape(Data,4,[]);
    Data_P = Data(1:2,:);
    Data_W = Data(3:4,:);

    for i=1:run
        pro = AttackMOKP('N',n,'maxFE',maxfe,'parameter',{Data_P,Data_W,0,0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end},sum(Data_P,2)');
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file,'origin');

elseif Pro == "PO"

    for i=1:run
        pro = AttackPO('N',n,'maxFE',maxfe,'parameter',{Data,0,0});
        Alg.Solve(pro);
        origin(i) = pro.CalMetric('HV',Alg.result{end});
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file,'origin');

elseif Pro == "OPF"
    
    for i=1:run
        pro = AttackCase118('N',n,'maxFE',maxfe,'parameter',{Data,0,0});
        Alg.Solve(pro);
        origin(i) = Min_value(Alg.result{end});
    end
    save(file,'origin');

elseif Pro == "CRT"

    for i=1:run
        pro = AttackCRT('N',n,'maxFE',maxfe,'parameter',{[],0,0});
        Alg.Solve(pro);
        origin(i) = HV(Alg.result{end},[1,2]);
        parsave(Alg.result{end}.objs, fileparts(save_foder), i, Pro, A_attacked);
    end
    save(file,'origin');
end
end


function parsave(objs, fa_file, i, Pro, A_attacked)
    file2             = sprintf('Val_%s_%s_%d.mat',Pro,A_attacked,i);
    file2             = fullfile(fa_file,file2);
    save(file2,"objs");
end