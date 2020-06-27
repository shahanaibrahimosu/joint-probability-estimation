function step_size=estimate_stepsize(A,lambda,marg,rho)
N= length(A);
F=size(A,2);
G=zeros(F,F);
s=zeros(N,1);
for n = 1:N
    for k=[1:n-1,n+1:N]
        for l=k+1:N
            G=G+(diag(lambda)*(A{l}'*A{l}).*(A{k}'*A{k})*diag(lambda));   
        end
    end
    s(n)=max(sum(G));
end
L=max(s);
step_size=rho/L;
end