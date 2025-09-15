function [Data, Encoding, Lower, Upper] = LoadOrgData(config)

    if config.Pro == "TSP"
        
        Datatemp   = load('Data_TSP-D200.mat');
        Data       = Datatemp.Data;

        Encoding   = ones(1,numel(Data));
        Lower      = -reshape(Data,1,[]);    
        Upper      = 1 - reshape(Data,1,[]);

    elseif config.Pro == "KP"

        load("P-100D.mat",'P');
        load("W-100D.mat",'W');
        Data       = [P,W];
        
        Encoding   = ones(1,numel(Data)); 
        Lower      = - reshape(Data,1,[]);     
        Upper      = 1 - reshape(Data,1,[]);

    elseif config.Pro == "MOTSP"

        load("MOTSP-D90.mat",'C');
        Data       = C;

        Encoding   = ones(1,2*10*(10-1)/2);
        index1     = find(tril(Data{1},-1)~=0);
        index2     = find(tril(Data{2},-1)~=0);
        temp1      = tril(Data{1},-1);
        temp2      = tril(Data{2},-1);
        Lower      = -[temp1(index1);temp2(index2)]';   
        Upper      = 1 - [temp1(index1);temp2(index2)]';


    elseif config.Pro == "MOKP"

        load("MOP-M2-D250.mat",'MOP');
        load("MOW-M2-D250.mat",'MOW');
        Data       = [MOP;MOW]; 
        Data       = reshape(Data,1,[]);

        Encoding   = ones(1,numel(Data));
        Lower      = - reshape(Data,1,[]);     
        Upper      = 1 - reshape(Data,1,[]);


    elseif config.Pro == "PO"

        load('Dataset_PO.mat','Dataset');
        Data       = Dataset.('data1000');
        Data       = Data(1:20,1:10);
        Datamin    = min(Data,[],'all');
        Datamax    = max(Data,[],'all');
        Datanorm   = (Data - Datamin)/(Datamax - Datamin);

        Encoding   = ones(1,numel(Datanorm));     
        Lower      = -reshape(Datanorm,1,[]);     
        Upper      = 1-reshape(Datanorm,1,[]);

    elseif config.Pro == "OPF"

        Data       = load('mpc.mat');
        Data       = Data.mpc.gencost(:,5:6);
        max_1      = max(Data(:,1));
        min_1      = min(Data(:,1));
        max_2      = max(Data(:,2));
        min_2      = min(Data(:,2));
        Datanorm(:,1) = (Data(:,1) - min_1)/(max_1 - min_1);
        Datanorm(:,2) = (Data(:,2) - min_2)/(max_2 - min_2);

        Encoding   = ones(1,numel(Data));     
        Lower      = -reshape(Datanorm,1,[]);     
        Upper      = 1-reshape(Datanorm,1,[]);

    elseif config.Pro == "CRT"

        load('Eyeball_L_PoiIndex_weight.mat');
        Eyeball_L  = Eyeball_L_PoiIndex_weight;
        nonzp      = find(Eyeball_L ~= 0);
        zp         = find(Eyeball_L == 0);
        num_to_select = 500 - numel(nonzp);
        szp        = randsample(zp,num_to_select);

        sze        = Eyeball_L(szp);

        sp         = [nonzp;szp];

        Data       = Eyeball_L(nonzp);

        Datanorm   = (Data - min(Data))/(max(Data)-min(Data));
        
        Encoding   = ones(1,numel(Data));     
        Lower      = -config.scalefactor*Datanorm';    
        Upper      = config.scalefactor*(1 - Datanorm');

        save(fullfile(fullfile(fileparts(fileparts(mfilename("fullpath"))),'PlatEMO 4.0'),'sp.mat'),'nonzp');
    else
        
        disp("please input Pro = TSP or KP or MOTSP or MOKP or PO or OPT or CRT");
        return;
        
    end
end

