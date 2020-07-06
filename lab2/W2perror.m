function res = W2perror(A,p,h)
N = 2/h+1;
U = reshape(A,[N N]);
dx = [1,0,1,1];
dy = [0,1,1,-1];
res = 0;
x = 2:N-1;
y = x;
for i = 1:4
    err = (U(x+dx(i),y+dy(i)) + U(x-dx(i),y-dy(i)) - 2*U(x,y))/((dx(i)^2+dy(i)^2)*h^2);
    res = res + sum(sum(abs(err).^p));
end
res = (res * h^2*3)^(1/p);
            
            