function U = WideStencil(U,F, Basis, g1, iter, tol)
M = size(U,1);
%% add the boundary constraint
eps  = 1e-8;
g.l  = @(y) g1.l(y).*(y>= eps).*(y<=1+eps);
g.r = @(y) g1.r(y).*(y>=-eps).*(y<=1-eps);
g.u =@(y) g1.u(y).*(y>=eps).*(y<=1+eps);
g.d = @(y) g1.d(y).*(y>=-eps).*(y<=1-eps);
LOOP = 10;

%Basis = {{[1,0],[0,1]},{[1,1],[-1,1]}};  % 9stencil
Res = reshape(MAFunction(U,Basis,g)-F, [M^2,1]);
for i = 1: iter
    Jac = @(v)(reshape(MAJacobi(U,Basis,g,reshape(v,[M,M])),[M^2,1]));
    [DU, flag] = gmres(Jac, Res, 1600);
    lambda =1;
    sigma = 0.6;
    innerloop = 0;
    if flag == 0
    DU = reshape(DU,[M,M]);
    disp(norm(DU,'fro'));
    DU = DU/(1+norm(DU,'fro'));
    Unew = U - lambda*DU;
    Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
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
    else
        Unew = (1+randn()/10)*U;
        Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
    end
%     if  innerloop > 2* LOOP 
%         Unew = U + rand()*F;%rand()/(M+1)^2*((1:M)'.^2+(1:M).^2);
%         Resnew = reshape(MAFunction(Unew,Basis,g)-F, [M^2,1]);
%     end

    
    disp(innerloop);
    U = Unew; 
    fprintf("Iter : %d: %e\n", i, norm(Resnew,'fro'));
    if (norm(Res,'fro')<tol*(1+norm(F,'fro'))); break; end
    Res = Resnew;
end

    
    


end

