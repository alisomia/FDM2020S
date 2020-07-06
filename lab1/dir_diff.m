function [du, rho] = dir_diff(u, i,j, g)
%U the array [i,j] the direction g the boudary value
   [M,~] = size(u);
   du = -u;
   rho = ones(M);
   xl = max(1-i,1); xr = min(M-i,M);
   yl = max(1-j,1); yr = min(M-j,M);
   % fill the center
    du(xl:xr,yl:yr) = du(xl:xr,yl:yr)+u(i+(xl:xr),j+(yl:yr));
    % fill right
    I = (xr+1:M)'; J = (1:M);
    du(I,J) = du(I,J) + g.r((J+(M+1-I).*(j/i))/(M+1));
    rho(I,J) = min(rho(I,J),repmat((M+1-I)./i,1,M));
    % fill up
    I = (1:M)'; J = (yr+1:M);
    du(I,J) = du(I,J) + g.u((I+(M+1-J).*(i/j))/(M+1));
    rho(I,J) = min(rho(I,J),repmat((M+1-J)./j, M,1));
    % fill left 
    I = (1:xl-1)'; J = (1:M);
    du(I,J) = du(I,J) + g.l((J - I.*(j/i))/(M+1));
    rho(I,J) = min(rho(I,J),repmat(I./abs(i),1,M));
    % fill down
    I = (1:M)'; J = (1:yl-1);
    du(I,J) = du(I,J) + g.d((I-J.*(i/j))/(M+1));
    rho(I,J) = min(rho(I,J),repmat(J./abs(j),M,1));
end
    
    
