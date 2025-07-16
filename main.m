function main(varargin)
    
    % Load parameters
    parser              = inputParser;

    addOptional(parser,"maxnumber",1);

    addOptional(parser,"Pro",'TSP',@ischar);
    addOptional(parser,"A_attacked",'GA',@ischar);
    addOptional(parser,"A_attack",'SparseAttack',@ischar);
    
    addOptional(parser,"maxfe",2500,@(x) x>0);
    addOptional(parser,"MAXFE1",3000,@(x) x>0);
    addOptional(parser,"n",30,@(x) x>0);
    addOptional(parser,"N1",30,@(x) x>0);

    addOptional(parser,"beta",0.05);
    addOptional(parser,"n_s",3);

    addOptional(parser,"t_s",0.1);
    
    addOptional(parser,"run",10);
    
    addOptional(parser,"scalefactor",1);

    parse(parser,varargin{:});
    
    config.maxnumber   = parser.Results.maxnumber;

    config.Pro         = parser.Results.Pro;
    config.A_attacked  = parser.Results.A_attacked;
    config.A_attack    = parser.Results.A_attack;

    config.maxfe       = parser.Results.maxfe;
    config.MAXFE       = parser.Results.MAXFE1;
    config.n           = parser.Results.n;
    config.N           = parser.Results.N1;

    config.beta        = parser.Results.beta;
    config.n_s         = parser.Results.n_s;

    config.t_s         = parser.Results.t_s;

    config.run         = parser.Results.run;
    
    config.scalefactor = parser.Results.scalefactor;

    mainfolder         = fileparts(mfilename("fullpath"));
    addpath(genpath(mainfolder));

    save_folder_pro    = fullfile(mainfolder,config.Pro);

    save_folder        = fullfile(save_folder_pro,sprintf("test%d",config.maxnumber));
    save_folder        = fullfile(save_folder,sprintf("%s_%s",config.A_attacked,config.A_attack));
    if ~exist(save_folder,'dir')
        mkdir(save_folder);
    end

    save(fullfile(save_folder,"config.mat"), '-struct', 'config');

    % load the raw data of the dataset
    [Data, Encoding, Lower, Upper] = LoadOrgData(config);

    origin            = Origin(config.Pro,Data,config.A_attacked,config.n,config.maxfe,config.run,save_folder,'save',0);
    origin            = mean(origin);

    Fun1              = @(x)mean(x~=0);
    Fun2              = Perturbation(config.Pro, Data, config.A_attacked, origin, config.n, config.maxfe);
    
    for fileno = 1:config.run
        [DEC_attack,OBJ_attack,~] = platemo('algorithm',{str2func(config.A_attack),config.beta,config.n_s,config.t_s},'N',config.N,'maxFE',config.MAXFE,'save',10,...
        'objFcn',{Fun1,Fun2},'encoding',Encoding,'lower',Lower,'upper',Upper);
        parsave(DEC_attack,OBJ_attack,config.A_attacked,config.A_attack,fileno, save_folder);
    end
    
end

function parsave(DEC_attack,OBJ_attack,AttackedAlgorithm,AttackAlgorithm,fileno,save_folder)
    file2             = sprintf('DEC_attack_%s_%s_%d.mat',AttackedAlgorithm,AttackAlgorithm,fileno);
    file2             = fullfile(save_folder,file2);
    save(file2,'DEC_attack');
    file3             = sprintf('OBJ_attack_%s_%s_%d.mat',AttackedAlgorithm,AttackAlgorithm,fileno);
    file3             = fullfile(save_folder,file3);
    save(file3,'OBJ_attack');
end