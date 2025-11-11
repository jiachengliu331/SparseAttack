function Obj2 = computeOriginObjSingle(Alg, Data, calperobj)
% COMPUTEORIGINOBJSINGLE - compute the objective for the best solution
% returned by an algorithm run (single-objective case).
%
% This is a renamed, clearer version of the old calOriginObj_Single.
% Keeping the behavior identical: take the first decision vector from
% Alg.result{end}.decs and evaluate it with calperobj.

    % Retrieve the decision matrix from the last recorded result
    dec = Alg.result{end}.decs;

    % Sanity check: ensure there is at least one decision vector
    if isempty(dec)
        error('computeOriginObjSingle:EmptyDecisions', 'Alg.result{end}.decs is empty.');
    end

    % Take the first (best) decision vector and compute its objective
    decmin = dec(1, :);
    Obj2 = calperobj(decmin, Data);

end
