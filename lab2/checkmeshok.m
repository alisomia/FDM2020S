function isok = checkmeshok(A,B)
ax = 1;
bx = find(B(1,:,1)==A(1,ax,1),1);
x = 1;
while 1
    ax = A(1,ax,3);
    bx = B(1,bx,3);
    if A(1,ax,1)~=B(1,bx,1)
        isok=0;
        break;
    end
    isok = 1;
    if ax==1
        return;
    end
end

