U = poissoninit(F,g);
V{1} = U;
global delta
for i = 2:8
    delta = 0.1^i*100;
    V{i} = WideStencil(V{i-1},F,Basis,g,100,1e-5);
end
norm(MAFunction(V{i},Basis,g)-F, 'fro')
U = V{i};