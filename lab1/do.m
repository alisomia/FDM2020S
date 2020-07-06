global delta
delta = 1e-6;
warning off
SIZE = [10,20,40,80];
[~,U0,F,g,Basis] = LoadFunction(1,1, 20);
WideStencil(poissoninit(F.*0+1,g),F,Basis,g,100,1e-5);
err = max(max(abs(U-U0)));
% for i = 4
%     for j = 1:3
%     [~,U0,F,g,Basis] = LoadFunction(3,j, SIZE(i));
%     [U, time] = solver(F,g,Basis);
%     err = max(max(abs(U-U0)));
%     RESULT(j, 2*i-1) = time;
%     RESULT(j,2*i) = err;
%     end
% end