function [U,F, G, X,Y,bd, elem, elemind, adj, count, bdadj] = getequation(u,f,h)
% Get the equation for Oliker--Prussner method.
%  
% INPUT
% ----------------------------------
% u [function (array, array) -> array]: the exact solution
% f [function (array, array) -> array]: det D^2u, for constructing RHS
% h [double]: mesh size
%
% OUTPUT
% ---------------------------------
% U [M*M double]: inital value
% F [M*M double]: RHS
% G [M*M double]: exact solution
% X [N double]: x-axis value of nodes
% Y [N double]: y-axis value of nodes
% bd [N logical]: 1 for boundary point and 0 for interior point                  
% elem [NT*3 double]: element array
%   - elem[i][1] elem[i][2] elem[i][3]: three nodes of element i,
%   counterclockwise. 
% elemind [NT*3 double]: element index array
%   - elemind[i][j]: index in the adjacent element of node elem[i][j] , means that
%   adj[elem[i][j]][elemind[i][j]][1] == i, used in fliping
% adj [N * L * 3 double]: adjacent linked list
%   - adj[i][k][1] the k-th element of adjacent elem of node i
%   - adj[i][k][2] the index of node in element adj[i].T[k][1], means that
%   elem[adj[i][k][1]][adj[i][k][2]] == i
%   - adj[i][k][3] the next element of adj[i].T[k][1], clockwise
% count [N double]: size for adjacent linked list of i
% bdadj [N*2 double]: neighborhood boundary point for all boundary point
%                               but corner one
%
%
% Linting@PKU
% 2020.06
X = -1:h:1;
Y = -1:h:1;
N = size(X,2);
[X,Y] = meshgrid(X,Y);
X = X(:);
Y = Y(:);
bd = (X.^2==1)|(Y.^2==1);
[elem, elemind, adj, count] = meshinit(X,Y,bd);

bdadj = zeros(size(bd,1),2);
bdadj(:,1) = (1:N^2)-N;
bdadj(:,2) = (1:N^2)+N;
bdadj(1:N, 1) = (1:N) - 1;
bdadj(1:N,2) = (1:N)+1;
bdadj(N^2-N+1:N^2,1) = (N^2-N+1:N^2) - 1;
bdadj(N^2-N+1:N^2,2) = (N^2-N+1:N^2) + 1;

[lambda,weight] = quadpts(4);
G = u(X,Y);
F = zeros(size(X));
for p = 1:size(lambda,1)
    Xp = lambda(p,1)*X(elem(:,1)) + lambda(p,2)*X(elem(:,2)) + lambda(p,3)*X(elem(:,3));
    Yp = lambda(p,1)*Y(elem(:,1)) + lambda(p,2)*Y(elem(:,2)) + lambda(p,3)*Y(elem(:,3));
    fp = f(Xp,Yp);
    for n = 1:size(elem,1)
    for j=1:3
        F(elem(n,j)) = F(elem(n,j)) + fp(n)*lambda(p,j)*h^2/2*weight(p);
    end
    end
end
U = (X.^2 + Y.^2 - 4)*4;
end
