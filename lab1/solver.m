% MA solver using Wide Stencil
% Input 
%--------------------------
% F[M*M double] RHS
% g[struct] BC, see LoadFunction
% Basis[Cell], Basis in wide stencil, see LoadFunction
%
% Output 
% -----------------------------
% U[M*M double] result
%
% delta can be modified, and used in max_delta and so on.

function [U, time] = solver(F, g, Basis)
tic
U = poissoninit(F,g);
V{1} = U;
global delta
for i = 2:9
    delta = 0.1^i*1000;
    V{i} = WideStencil(V{i-1},F,Basis,g,100,1e-5);
end
norm(MAFunction(V{i},Basis,g)-F, 'fro')
U = V{i};
time = toc; 
end