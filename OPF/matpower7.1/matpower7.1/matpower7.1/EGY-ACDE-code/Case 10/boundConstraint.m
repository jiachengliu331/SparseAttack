function vi = boundConstraint (vi, pop, X_max, X_min)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound

[NP, ~] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(X_min, NP, 1);
pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;

%% check the upper bound
xu = repmat(X_max, NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;