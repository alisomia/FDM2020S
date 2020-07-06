[u,f] = loadfunction(3);
h = 1/4;
[U,F, G, X,Y,bd, elem, elemind, adj, count, bdadj] = getequation(u,f,h);
[U, X,Y, bd, elem,elemind, adj, count] = OPinit(U,F, G, X,Y,bd, elem, elemind, adj, count,1000,bdadj, 1e-6);
tic
[U, X,Y, bd, elem,elemind, adj, count]  = OPsolve(U,F, G,X,Y,bd,elem, elemind, adj, count, 30000,1e-8);
toc

% [u,f] = loadfunction(3);
% H = [1,1/2,1/4,1/8];
% for i = 1:4
% h = H(i);
% [U,F, G, X,Y,bd, elem, elemind, adj, count, bdadj] = getequation(u,f,h);
% tic
% [U, X,Y, bd, elem,elemind, adj, count] = OPinit(U,F, G, X,Y,bd, elem, elemind, adj, count,100000,bdadj, 1e-6);
% liftbdtime = toc;
% tic
% [U, X,Y, bd, elem,elemind, adj, count]  = OPsolve(U,F, G,X,Y,bd,elem, elemind, adj, count, 100000,1e-6);
% RESULT(i,1)=toc;
% RESULT(i,2) = norm(G-U,'inf');
% RESULT(i,4)= W2perror(G-U,2,h);
% RESULT(i,6) = W2perror(G-U,3,h);
% end
% for i = 2:4
% for j = [2,4,6]
% RESULT(i,j+1) = log2(RESULT(i-1,j)/RESULT(i,j));
% end
% end