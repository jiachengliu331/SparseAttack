classdef UserProblem < PROBLEM
%UserProblem - User defined problem.
%
%   This problem provides a general framework, whose details can be defined
%   by the inputs of the constructor.
%
% All the acceptable properties:
%   encoding  	<string>            encoding scheme of each decision variable (1.real 2.integer 3.label 4.binary 5.permutation)
%   lower    	<vector>            lower bound of each decision variable
%   upper      	<vector>            upper bound of each decision variable
%   initFcn     <function handle>   function for initializing solutions
%   evalFcn     <function handle>   function for evaluating solutions
%   decFcn      <function handle>   function for repairing invalid solutions
%   objFcn     	<function handle>   objective functions
%   conFcn     	<function handle>   constraint functions
%   objGradFcn  <function handle>   function for calculating the gradients of objectives
%   objConFcn   <function handle>   function for calculating the gradients of constraints
%   parameter   <any>               dataset of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        initFcn    = {};        % Function for initializing solutions
        evalFcn    = {};        % Function for evaluating solutions
        decFcn     = {};    	% Function for repairing invalid solutions
        objFcn     = {};     	% Objective functions
        conFcn     = {};     	% Constraint functions
        objGradFcn = {};        % Function for calculating the gradients of objectives
        conGradFcn = {};        % Function for calculating the gradients of constraints
    end
    methods
        %% Constructor
        function obj = UserProblem(varargin)
            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','M','D','maxFE','maxRuntime','encoding','lower','upper','initFcn','evalFcn','decFcn','objFcn','conFcn','objGradFcn','conGradFcn','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            if isempty(obj.D)
                obj.encoding = Str2Fcn(obj.encoding,1,[],'encoding scheme');
                obj.D = length(obj.encoding);
            else
                obj.encoding = Str2Fcn(obj.encoding,1,[],'encoding scheme',obj.D);
            end
            obj.lower      = Str2Fcn(obj.lower,1,[],'lower bound',obj.D);
            obj.upper      = Str2Fcn(obj.upper,1,[],'upper bound',obj.D);
            obj.parameter  = Str2Fcn(obj.parameter,1,[],'dataset');
            obj.initFcn    = Str2Fcn(obj.initFcn,2,~isempty(obj.parameter),'initialization function');
            obj.evalFcn    = Str2Fcn(obj.evalFcn,3,~isempty(obj.parameter),'evaluation function');
            obj.decFcn     = Str2Fcn(obj.decFcn,3,~isempty(obj.parameter),'repair function');
            obj.objFcn     = Strs2Fcns(obj.objFcn,3,~isempty(obj.parameter),'objective function f');
            obj.conFcn     = Strs2Fcns(obj.conFcn,3,~isempty(obj.parameter),'constraint function g');
            obj.objGradFcn = Strs2Fcns(obj.objGradFcn,3,~isempty(obj.parameter),'gradient of objective fg');
            obj.conGradFcn = Strs2Fcns(obj.conGradFcn,3,~isempty(obj.parameter),'gradient of constraint gg');
            if isempty(obj.evalFcn)
                obj.M = length(obj.objFcn);
            else
               [~,f,~] = CallFcn(obj.evalFcn,1,obj.parameter);
               obj.M   = length(f);
            end
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            if ~isempty(obj.initFcn)
                Population = obj.Evaluation(CallFcn(obj.initFcn,N,obj.parameter,'initialization function',[N obj.D]));
            else
                Population = Initialization@PROBLEM(obj,N);
            end
            obj.optimum = max(Population.objs,[],1);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            if ~isempty(obj.evalFcn)
                for i = 1 : size(varargin{1},1)
                    [PopDec(i,:),PopObj(i,:),PopCon(i,:)] = CallFcn(obj.evalFcn,varargin{1}(i,:),obj.parameter,'evaluation function',[1 obj.D],[1 obj.M]);
                end
                Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
                obj.FE     = obj.FE + length(Population);
            else
                Population = Evaluation@PROBLEM(obj,varargin{:});
            end
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            if ~isempty(obj.decFcn)
                for i = 1 : size(PopDec,1)
                    PopDec(i,:) = CallFcn(obj.decFcn,PopDec(i,:),obj.parameter,'repair function',[1 obj.D]);
                end
            else
                PopDec = CalDec@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            if ~isempty(obj.objFcn)
                PopObj = zeros(size(PopDec,1),length(obj.objFcn));
                for i = 1 : size(PopDec,1)
                    for j = 1 : length(obj.objFcn)
%                         PopObj(i,j) = CallFcn(obj.objFcn{j},PopDec(i,:),obj.parameter,sprintf('objective function f%d',j),[1 1]);
                          if j==length(obj.objFcn)
                              [val1,val2] = CallFcn(obj.objFcn{j},PopDec(i,:),obj.parameter,sprintf('objective function f%d',j),[1 1]);
                              PopObj(i,j) = val1;
                              obj.val2{i,1}   = val2;
                          else
                              PopObj(i,j) = CallFcn(obj.objFcn{j},PopDec(i,:),obj.parameter,sprintf('objective function f%d',j),[1 1]);
                          end
                    end
                end
            else
                PopObj = CalObj@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            if ~isempty(obj.conFcn)
                PopCon = zeros(size(PopDec,1),length(obj.conFcn));
                for i = 1 : size(PopDec,1)
                    for j = 1 : length(obj.conFcn)
                        PopCon(i,j) = CallFcn(obj.conFcn{j},PopDec(i,:),obj.parameter,sprintf('constraint function g%d',j),[1 1]);
                    end
                end
            else
                PopCon = CalCon@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate the gradients of objectives
        function ObjGrad = CalObjGrad(obj,Dec)
            if ~isempty(obj.objGradFcn)
                ObjGrad = zeros(length(obj.objGradFcn),obj.D);
                for i = 1 : length(obj.objGradFcn)
                    ObjGrad(i,:) = CallFcn(obj.objGradFcn{i},Dec,obj.parameter,sprintf('gradient of objective fg%d',i),[1 obj.D]);
                end
            else
                ObjGrad = CalObjGrad@PROBLEM(obj,Dec);
            end
        end
        %% Calculate the gradients of constraints
        function ConGrad = CalConGrad(obj,Dec)
            if ~isempty(obj.conGradFcn)
                ConGrad = zeros(length(obj.conGradFcn),obj.D);
                for i = 1 : length(obj.conGradFcn)
                    ConGrad(i,:) = CallFcn(obj.conGradFcn{i},Dec,obj.parameter,sprintf('gradient of constraint gg%d',i),[1 obj.D]);
                end
            else
                ConGrad = CalConGrad@PROBLEM(obj,Dec);
            end
        end
    end
end

function var = Str2Fcn(var,type,useData,name,D)
% Convert a string into a variable or function and check its validity

    if ischar(var)
        try
            if ~isempty(regexp(var,'^<.+>$','once'))
                switch type
                    case 1
                        var = load(var(2:end-1));
                    otherwise
                        [folder,file] = fileparts(var(2:end-1));
                        addpath(folder);
                        var = str2func(file);
                end
            else
                switch type
                    case 1
                        var = str2num(var);
                    case 2
                        if useData
                            var = str2func(['@(N,data)',var]);
                        else
                            var = str2func(['@(N)',var]);
                        end
                    case 3
                        if useData
                            var = str2func(['@(x,data)',var]);
                        else
                            var = str2func(['@(x)',var]);
                        end
                end
            end
        catch err
            err = addCause(err,MException('','Fail to define the %s',name));
            rethrow(err);
        end
    end
    if type == 1 && nargin > 4
        if isscalar(var)
            var = repmat(var,1,D);
        else
            assert(ismatrix(var)&&all(size(var)==[1,D]),'the %s should be a scalar or a 1*%d vector, while its current size is %d*%d.',name,D,size(var,1),size(var,2));
        end
    end
end

function Var = Strs2Fcns(Var,type,useData,name)
% Convert multiple strings into functions

    if ~iscell(Var)
        Var = {Var};
    end
    Var(cellfun(@isempty,Var)) = [];
    for i = 1 : length(Var)
        Var{i} = Str2Fcn(Var{i},type,useData,[name,num2str(i)]);
    end
end

function varargout = CallFcn(func,input,parameter,name,varargin)
% Call a function and check the validity of its output

    try
        if isempty(parameter)
            [varargout{1:nargout}] = func(input);
        else
            [varargout{1:nargout}] = func(input,parameter);
        end
        for i = 1 : min(length(varargout),length(varargin))
            assert(ismatrix(varargout{i})&&all(size(varargout{i})==varargin{i}),'the size of its output #%d should be %d*%d, while its current size is %d*%d.',i,varargin{i}(1),varargin{i}(2),size(varargout{i},1),size(varargout{i},2));
        end
    catch err
        err = addCause(err,MException('','The %s is invalid',name));
        rethrow(err);
    end
end