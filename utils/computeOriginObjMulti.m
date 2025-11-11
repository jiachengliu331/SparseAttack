function [HV2,obj_f1] = computeOriginObjMulti(Alg, Data, calperobj1, points)
% COMPUTEORIGINOBJMULTI - compute multi-objective performance (HV) of
% the final population returned by an algorithm run.
%
% Syntax:
%   [HV2, obj_f1] = computeOriginObjMulti(Alg, Data, calperobj1, points)
%
% Inputs:
%   Alg        - algorithm object/struct with Alg.result containing final population
%   Data       - original problem data passed to evaluator
%   calperobj1 - function handle: @(PopDec, Data) -> objective matrix for population
%   points     - reference points used by HV test (depends on problem)
%
% Outputs:
%   HV2        - hypervolume value computed from the evaluated population
%   obj_f1     - objective matrix for the evaluated population

    % decp_f1: obtain the population decision matrix used for objective 1
    decp_f1 = findf1(Alg.result{end});

    % compute objective values for the population using provided evaluator
    obj_f1 = calperobj1(decp_f1, Data);

    % compute hypervolume (or other multi-objective metric provided by HVtest)
    HV2 = HVtest(obj_f1, points);

end
