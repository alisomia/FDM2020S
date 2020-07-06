% Using Wide stencil method to solve monge ampere equation via Damped
% Newton method. 
% 
% input
% ----------------------------------------
% U[M*M double]: initial value
% F[M*M double]: RHS
% Basis[cell]: see LoadFunction
% g1[struct]: BC, see LoadFunction
% iter[int] iteration time 
% tol[double]: tolerance
function U = WideStencil(U,F, Basis, g1, iter, tol)
M = size(U,1);
%% add the boundary constraint to avoid repeating
eps  = 1e-8;
g.l  = @(y) g1.l(y).*(y>= eps).*(y<=1+eps);
g.r = @(y) g1.r(y).*(y>=-eps).*(y<=1-eps);
g.u =@(y) g1.u(y).*(y>=eps).*(y<=1+eps);
g.d = @(y) g1.d(y).*(y>=-eps).*(y<=1-eps);
LOOP = 10;  % inexact linesearch time


Res = reshape(MAFunction(U,Basis,g)-F, [M^2,1]); %Compute the residue
for i = 1: iter
    Jac = @(v)(reshape(MAJacobi(U,Basis,g,reshape(v,[M,M])),[M^2,1]));
    [DU, ~] = gmres(Jac, Res, 1600);
    lambda =1;
    sigma = 0.6;
    innerloop = 0;
    DU = reshape(DU,[M,M]);
    %sum(sum(Jac(DU).*Res))
    %disp(norm(DU,'fro'));
    DU = DU/(1+norm(DU,'fro'));
    Unew = U - lambda*DU;
    Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
    % START LINE SEARCH
    while norm(Resnew,'fro') > norm(Res,'fro')
        lambda = lambda*sigma;
        Unew = U - lambda*DU;
        Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
        innerloop = innerloop + 1;
       if innerloop >LOOP; break; end
    end
    if innerloop>LOOP
         lambda = 1;
         Unew = U + lambda*DU;
        Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
       
        while norm(Resnew,'fro') > norm(Res,'fro')
            lambda = lambda*sigma;
            Unew = U + lambda*DU;
            Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
            innerloop = innerloop + 1;
            if innerloop >2*LOOP; break; end
        end
    end
    % END LINE SEARCH
    
%     if abs(norm(Resnew, 'fro') - norm(Res, 'fro')) < 1e-4
%         Unew = U + rand()/(M+1)^2*((1:M)'.^2+(1:M).^2);
%         Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
%     end
%     if  innerloop > 2* LOOP 
%         Unew = U + rand()*F;%rand()/(M+1)^2*((1:M)'.^2+(1:M).^2);
%         Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
%     end

    U = Unew; 
    fprintf("Iter : %d: %e\n", i, norm(Resnew,'fro'));
    Res = Resnew;
    if (norm(Res,'fro')<tol*(1+norm(F,'fro'))); break; end
    
end

    
    


end

