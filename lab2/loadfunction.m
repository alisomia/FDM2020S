function [u,f] = loadfunction(id)
% load equation for given id
% u is the exact solution, f = det D^2u
switch id
    case 0 
        u = @(x,y) x.^2 + y.^2 +1;
        f = @(x,y) 4+x.*0;
    case 1
        u = @(x,y) exp((x.^2+y.^2)/2);
        f = @(x,y) exp((x.^2+y.^2)).*(1+x.^2+y.^2);
    case 2
        u = @(x,y) 2*max(sqrt(x.^2+y.^2)-0.5,0).^2 + 2*(x.^2+y.^2);
        f = @(x,y) 16*(x.^2+y.^2<1/4) + (64-16./sqrt(x.^2+y.^2)).*(x.^2+y.^2>=1/4);
    case 3
        u1 = @(x,y) (x.^4 + 1.5*(x+1e-8).^(-2).*y.^2).*(abs(y)<=abs(x.^3)) + (0.5*x.^2.*y.^(2/3)+2*y.^(4/3)).*(abs(y)>abs(x.^3));
        f1 = @(x,y) (36-9*x.^(-6).*y.^2).*(abs(y)<=abs(x.^3)) + (8/9-5/9*x.^2.*y.^(-2/3)).*(abs(y)>abs(x.^3));
        u = @(x,y) u1(abs(x), abs(y));
        f = @(x,y) f1(abs(x),abs(y));
end
end