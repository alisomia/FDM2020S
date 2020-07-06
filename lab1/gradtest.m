eps = 1e-6;
S = randn(9)*0.01;
JS = MAJacobi(U, Basis, g, S);
JJS = (MAFunction(U+eps*S, Basis, g) - MAFunction(U-eps*S,Basis,g))/eps/2;
norm(JS-JJS,'fro')