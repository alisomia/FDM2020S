function [U, X,Y, bd, elem,elemind, adj, count] = OPinit(U,~, G, X,Y,bd,elem, elemind, adj, count, maxiter, bdadj,tol)
% lift the boundary
% 
% Linting@PKU
% 2020.06
next = [2,3,1];
inf = 1e10;
eps = 1e-6;

for iter = 1:maxiter
    for i = 1:size(count)
        if ~bd(i); continue; end  % skip interior point
        if G(i) < U(i)+eps; continue; end % skip points lifted already
        while 1
            %%% begin computing threshold
            P = zeros(count(i)+1,2);
            Q = P; N = P;
            th_adj = ones(10,1)*inf; % stores the threshold, init with inf
            t2 = find(adj(i,:,1),1);
            cnt = 1;
            [P(cnt,:),Q(cnt,:),N(cnt,:)] = getdis(i,t2);
            
            % traversal the linked list
            for cnt = 1:count(i)
                % compute the threshold for each edge adj to a
                t1 = t2;
                t2 = adj(i,t2,3);
                j = elem(adj(i,t1,1),next(adj(i,t1,2)));
                [P(cnt+1,:),Q(cnt+1,:),N(cnt+1,:)] = getdis(i,t2);
                if bd(j) && ((X(i)==X(j))||(Y(i)==Y(j)))  % skip the virtual link
                    continue;
                end
                th_adj(t1) = (N(cnt,:)*P(cnt,:)'-N(cnt,:)*P(cnt+1,:)')/(N(cnt,:)*Q(cnt,:)' - N(cnt,:)*Q(cnt+1,:)');
                if th_adj(t1)<U(i); th_adj(t1) = inf; end
            end
            [th0_adj, t0] = min(th_adj);
            %%% end computing threshold
            
            % compute boundary threshold
            if bd(i)==1 && ~((X(i)^2==1)&&(Y(i)^2==1)) % judge if it is a corner point
                th0_line = 0.5*(U(bdadj(i,1))+U(bdadj(i,2)));
            else
                th0_line = inf;
            end
            
            th0 = min([th0_adj, th0_line]);
            
            if th0>G(i)
                U(i) = G(i);
                break
            end
            
            if th0>th0_line-eps
                U(i) = th0_line;
                break;
            else
                
                %%% begin edge flipping
                T0 = adj(i,t0,1); t1 = adj(i,t0,3); T1 = adj(i,t1,1);
                j = elem(T0, next(adj(i,t0,2)));
                k = elem(T0, next(next(adj(i,t0,2))));
                l = elem(T1, next(adj(i,t1,2)));
                
                % change the elem
                elem(T0,:) = [i,l,k];
                elem(T1,:) = [j,k,l];
                
                
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
                
                %%% end edge fliping
                if U(i) < th0+ eps
                    continue
                else
                    U(i) = th0;
                    break;
                end
            end
        end
        
    end
    res = (G-U)'*bd;
    %fprintf("residual=%f\n", res);
    if res<tol; return; end
end

    function [pp,qq,nn] = getdis(i,t)
        
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
        pp = [AAinv11*U(elem(T,b)) + AAinv12 * U(elem(T,c)),...
            AAinv21*U(elem(T,b)) + AAinv22 * U(elem(T,c))];
        qq = [AAinv11 + AAinv12; AAinv21 + AAinv22];
        nn = [AA12, -AA11];
    end

end