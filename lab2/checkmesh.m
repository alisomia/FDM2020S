function s = checkmesh(X,Y,bd,elem,i)
[~,~,adj,~] = meshinit(X,Y,bd,elem);
s = adj(i,:,:);