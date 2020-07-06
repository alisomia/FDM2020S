% Compute the Monge Ampere Discrete Variational Formula
function [F] = MAFunction(U, Basis, g)
NBasis = size(Basis,2);
MA = cell(1,NBasis);
for i = 1:NBasis
    A =central_diff(U,Basis{i}{1},g);
    B = central_diff(U,Basis{i}{2},g);
    MA{i} = max_delta(A,0).*max_delta(B,0);
end
F = min_delta(MA{1}, MA{2});
for i=3:NBasis
    F = min_delta(F, MA{i});
end
end

