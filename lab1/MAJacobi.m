% Compute the Directional Derivative of Monge Ampere Discrete Variational Formula
% The Jacobian can be used by Jac = @(v) MAJacobi(U, Basis, g, v)
% when Krylov Subspace Method is called.
function [G] = MAJacobi(U, Basis, g, V)
NBasis = size(Basis,2);

% zero bdry
gg.d = @(x) 0;
gg.u = @(x) 0;
gg.l = @(x) 0;
gg.r = @(x) 0;

MA = cell(1,NBasis);
MAJ = cell(1,NBasis);
for i = 1:NBasis
    A =central_diff(U,Basis{i}{1},g);
    B = central_diff(U,Basis{i}{2},g);
    GA = central_diff(V, Basis{i}{1}, gg);
    GB = central_diff(V, Basis{i}{2}, gg);
    MA{i} = max_delta(A,0).*max_delta(B,0);
    MAJ{i} = max_delta_grad(A,0).*GA.*max_delta(B,0) + max_delta_grad(B,0).*GB.*max_delta(A,0);
end
F = min_delta(MA{1}, MA{2});
G = min_delta_grad(MA{1}, MA{2}).*MAJ{1} + min_delta_grad(MA{2},MA{1}).*MAJ{2};
for i=3:NBasis
    G = min_delta_grad(F,MA{i}).*G + min_delta_grad(MA{i}, F).*MAJ{i};
    F = min_delta(F, MA{i});
end
end
