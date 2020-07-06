function [U, X,Y, bd, elem,elemind, adj, count] = OPsolve(U,F,G,X,Y,bd,elem, elemind, adj, count, maxiter, tol)
% solve oliker--prussner method by Perron's iteration
%
% Linting@PKU
% 202.06
next = [2,3,1];
inf = 1e10;
eps = 1e-6;
loss = zeros(maxiter:1);
for iter = 1:maxiter
    Uold = U;
    for i = 1:size(count,1)
        if bd(i); continue; end % skip boundary point
        while 1
            %%% begin computing threshold
            P = zeros(count(i)+1,2);
            Q = P; N = P;
            th_adj = ones(size(adj,2)+1,1)*inf; % stores the threshold
            t2 = find(adj(i,:,1),1);
            cnt = 1;
            [P(cnt,1),P(cnt,2),Q(cnt,1),Q(cnt,2),N(cnt,1),N(cnt,2)] = getdis(i,t2);
            
            
            for cnt = 1:count(i)
                % compute the threshold for each edge adj to a
                t1 = t2;
                t2 = adj(i,t2,3);
                [P(cnt+1,1),P(cnt+1,2),Q(cnt+1,1),Q(cnt+1,2),N(cnt+1,1),N(cnt+1,2)] = getdis(i,t2);
                th_adj(t1) = (N(cnt,1)*P(cnt,1) + N(cnt,2)*P(cnt,2) - N(cnt,1)*P(cnt+1,1) - N(cnt,2)*P(cnt+1,2))/....
                    (N(cnt,1)*Q(cnt,1) + N(cnt,2)*Q(cnt,2) - N(cnt,1)*Q(cnt+1,1) - N(cnt,2)*Q(cnt+1,2));
                if th_adj(t1)<U(i)-eps; th_adj(t1) = inf; end
            end
            
            
            [th0, t0] = min(th_adj);
            th0 = th0 - eps/100;
            
            %%% end threshold
            
            if subgrad(P,Q,th0)>F(i)+eps
                
                % edge flipping
                T0 = adj(i,t0,1); t1 = adj(i,t0,3); T1 = adj(i,t1,1);
                j = elem(T0, next(adj(i,t0,2)));
                k = elem(T0, next(next(adj(i,t0,2))));
                l = elem(T1, next(adj(i,t1,2)));
                
                % change elem
                elem(T0,:) = [i,l,k];
                elem(T1,:) = [j,k,l];
                
                
                % change adj
                ind_i_T0 = t0;
                ind_i_T1 = t1;
                ind_j_T1 = elemind(T1,next(next(adj(i,ind_i_T1,2))));
                ind_j_T0 = adj(j,ind_j_T1,3);
                ind_k_T0 = elemind(T0,next(next(adj(i,ind_i_T0,2))));
                ind_k_newT0 = find(adj(k,:,1)==0,1);
                
                ind_l_T1 = elemind(T1,next(adj(i,ind_i_T1,2)));
                ind_l_newT1 = find(adj(l,:,1)==0,1);
                
                adj(i,ind_i_T1,1) = 0;
                adj(i,ind_i_T0,3) = adj(i,ind_i_T1,3);
                adj(i,ind_i_T0,2) = 1;
                count(i) = count(i) - 1;
                
                adj(j,ind_j_T0,1) = 0;
                adj(j,ind_j_T1,3) = adj(j,ind_j_T0,3);
                adj(j,ind_j_T1,2) = 1;
                count(j) = count(j) - 1;
                
                adj(k,ind_k_newT0,3) = adj(k,ind_k_T0,3);
                adj(k,ind_k_T0,3) = ind_k_newT0;
                adj(k,ind_k_T0,1) = T1;
                adj(k,ind_k_T0,2) = 2;
                adj(k,ind_k_newT0,1) = T0;
                adj(k,ind_k_newT0,2) = 3;
                count(k) = count(k) + 1;
                
                adj(l,ind_l_newT1,3) = adj(l,ind_l_T1,3);
                adj(l,ind_l_T1,3) = ind_l_newT1;
                adj(l,ind_l_T1,1) = T0;
                adj(l,ind_l_T1,2) = 2;
                adj(l,ind_l_newT1,1) =T1;
                adj(l,ind_l_newT1,2) = 3;
                count(l) = count(l)+1;
                
                %change elemind
                elemind(T0,1) = ind_i_T0;
                elemind(T0,2) = ind_l_T1;
                elemind(T0,3) = ind_k_newT0;
                elemind(T1,1) = ind_j_T1;
                elemind(T1,2) = ind_k_T0;
                elemind(T1,3) = ind_l_newT1;
                disp('kao');
                U(i) = th0;
                
                continue;
            else

                ll = U(i); rr = th0 + eps/100;
                C = subgrad(P,Q,ll)/(ll-rr)^2;
                U(i) = rr - sqrt(F(i)/C);
                break;
            end
        end
 
    end

    fprintf("iter = %d, loss = %e, practical loss = %e\n", iter, norm(U-G,'inf'), norm(U-Uold,'inf'));
    %loss(iter) = norm(U-G,'inf');
    if norm(U-Uold,'inf')<tol; break; end
end
%plot(1:iter,log(loss(1:iter)));
disp(iter)
    function [pp1,pp2,qq1,qq2,nn1,nn2] = getdis(i,t)
        a = adj(i,t,2); b = next(a); c = next(b);
        T = adj(i,t,1);
%         AA = [X(elem(T,b)) - X(elem(T,a)),...
%             Y(elem(T,b)) - Y(elem(T,a));...
%             X(elem(T,c)) - X(elem(T,a)),...
%             Y(elem(T,c)) - Y(elem(T,a))];
%         pp = AA\[U(elem(T,b)); U(elem(T,c))];
%         qq = AA\[1;1];
%         nn = AA(1,:)';
%         nn = [nn(2),-nn(1)];
        AA11 = X(elem(T,b)) - X(elem(T,a));
        AA12 = Y(elem(T,b)) - Y(elem(T,a));
        AA21 = X(elem(T,c)) - X(elem(T,a));
        AA22 = Y(elem(T,c)) - Y(elem(T,a));
        detAA = AA11*AA22-AA12*AA21;
        AAinv11 = AA22/detAA;
        AAinv22 = AA11/detAA;
        AAinv12 = -AA12/detAA;
        AAinv21 = -AA21/detAA;
        pp1 = AAinv11*U(elem(T,b)) + AAinv12 * U(elem(T,c));
        pp2=   AAinv21*U(elem(T,b)) + AAinv22 * U(elem(T,c));
        qq1 = AAinv11 + AAinv12;
        qq2 = AAinv21 + AAinv22;
        nn1 = AA12;
        nn2 = -AA11;
    end

    function fi = subgrad(P,Q,ui)
        grad = P-ui*Q;
        fi = myarea(grad(:,1),grad(:,2));
    end
    function A = myarea(XX,YY)
        ii=2;
        A = abs((XX(ii)-XX(1))*(YY(ii+1)-YY(1))-(XX(ii+1)-XX(1))*(YY(ii)-YY(1)));
        for ii=3:size(XX,1)-1
            A = A+abs((XX(ii)-XX(1))*(YY(ii+1)-YY(1))-(XX(ii+1)-XX(1))*(YY(ii)-YY(1)));
        end
        A = A/2;
    end
end





