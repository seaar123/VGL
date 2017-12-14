function [a,s,Yhat]=NET(X,W,m_in,N,n_op,knd)
n=N+n_op+1;% no. of neurons ( include inputs )
a=zeros(1,n);% +1 for bais
a(1)=1;%#2 bais
a(2:m_in+1)=X;% #3
s=zeros(1,n);
for j=m_in+1+1:n
    s(j)=0;
    for m=1:j-1
        s(j)=s(j)+W(m,j)*a(m); % #5
    end
    [a(j),D_outs]=Activity_functions(s(j),knd);% #6
end
Yhat(1:n_op)=a(N+1+1:n);% #8
Yhat=Yhat';
end




