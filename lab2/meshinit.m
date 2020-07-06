% initialize the mesh with delaunay triangulation
%
%
% INPUT
% --------------------------------------------------
% X [N double]: x-axis value of nodes
% Y [N double]: y-axis value of nodes
% bd [N logical]: 1 for boundary point and 0 for interior point    
%
% OUTPUT
% ---------------------------------            
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
% Linting@PKU
% 2020.06
function [elem, elemind, adj, count] = meshinit(X,Y, ~, elem)
n = size(X,1);
% get delaunay triangultion 
if nargin<4; elem = delaunay(X,Y); end

% zero initialization
elemind = zeros(size(elem));
adj = zeros(n,40,3);
count = zeros(n,1);

% get adjacent array
for i = 1: size(elem,1)
    for j = 1:3
        k = elem(i,j);
        count(k) = count(k)+1;
        adj(k,count(k),1) = i;
        adj(k,count(k),2) = j;
        elemind(i,j) = count(k);
    end
end
next = [2,3,1];
% sort the adjacent array
for i = 1:size(count,1)
    cnt = 1:count(i);
    adjnode = elem(sub2ind(size(elem), adj(i,cnt,1), next(adj(i,cnt,2))));
    Xlocal = X(adjnode) - X(i);
    Ylocal = Y(adjnode) - Y(i);
    arg = atan2(Ylocal, Xlocal);
    [~,IX] = sort(arg, 'descend');
    adj(i,IX,3) = circshift(IX,-1);
end
    
end
