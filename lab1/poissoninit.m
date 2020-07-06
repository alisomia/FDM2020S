% get a proper initial value py solving poisson
% input
% ------------------------------
% F[M*M double] RHS
% g[struct] boundary condition, see LoadFunction
%
% output 
%----------------------------------
% U the result by solving lapalce u = 2sqrt(f), with BC g.
% U : [M*M double]
function U = poissoninit(F, g)
    M = size(F,1);
    F = 2*sqrt(F);
        X = (1:M)/(M+1);
    F(1,:) = F(1,:) - g.l(X)*(M+1)^2;
    F(M,:) = F(M,:) - g.r(X)*(M+1)^2;
    F(:,1) = F(:,1)  - g.d(X')*(M+1)^2;
    F(:,M) = F(:,M) - g.u(X')*(M+1)^2;
    %U = gmres(@(u) poi(u,g,M), ones(M^2,1), 1600,[], [], [],[],[] );
    U = minres(@(u) poi(u,g,M), reshape(F,[M^2,1]), [],1600 );
    U = reshape(U,[M,M]);
end
 %FF = [g.l(X); U(1:M-1,:)] + [U(2:M,:); g.r(X)] + [g.d(X'), U(:,1:M-1)] + [U(:,2:M), g.u(X')] - 4*U;

function FF = poi(U, g, M)
    U = reshape(U, [M,M]);
    %FF = [g.l(X'), U(:,1:M-1)] + [U(:,2:M), g.r(X')] + [g.u(X); U(1:M-1,:)] + [U(2:M,:); g.d(X)] - 4*U;
    FF = [zeros(1,M); U(1:M-1,:)] + [U(2:M,:); zeros(1,M)] + [zeros(M,1), U(:,1:M-1)] + [U(:,2:M), zeros(M,1)] - 4*U;
    FF = reshape(FF, [M^2,1])*(M+1)^2;
end


