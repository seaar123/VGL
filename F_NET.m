function [dy_dw,dy_dx]=F_NET(W,a,s,F_Yhat,m_in,N,n_op,knd)
n=N+n_op+1;% no. of neurons ( incluse inputs )
% F_W=zeros(n+N,n+N);
F_s=zeros(1,n);
% #11a
F_a=zeros(1,n);
%Getting dy/dw
for i=1:n_op
    F_a(i+N+1)=F_Yhat(i);
end
for j=n:-1:m_in+1+1 %1
    for m=j+1:n
        F_a(j)=F_a(j)+W(j,m)*F_s(m);% #11b
    end
    [outs,D_outs]=Activity_functions(s(j),knd);
    F_s(j)=F_a(j)*D_outs;% #12
    for m=1:j-1
        dy_dw(m,j)=F_s(j)*a(m);% #13
    end
end
%Getting dy/dx
for k=n_op:-1:1
    for i=1:n_op
        F_a(i+N+1)=F_Yhat(i);
    end
    F_a(k+N+1)=0;
    for i=m_in+1+1:n
        for j=1:i-1
            W(j,k+N+1)=0;
        end
    end
    for j=n:-1:1%m_in+1+1
        for m=j+1:n
            F_a(j)=F_a(j)+W(j,m)*F_s(m);
        end
        [outs,D_outs]=Activity_functions(s(j),knd);
        F_s(j)=F_a(j)*D_outs;
    end
    for j=2:m_in+1
        dy_dx(n_op-k+1,j-1)=F_a(j);
    end
end
    
    
    
