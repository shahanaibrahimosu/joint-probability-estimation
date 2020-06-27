function grad_array = computefullGradient(A,lambda,marg,J,J_order,Y)
N=size(A,1);
F=size(lambda,1);
grad_array = cell(N+1,1);
marg_array=cell2mat(marg);
for n=1:N
    [m,u]=size(A{n});
    grad_array{n}=zeros(m*u,J(n));
    for j=1:J(n)
        sel=marg_array(J_order{n}(j),:);
        ind=sel(sel~=n);
        X=Y{marg_array(:,1)==sel(1) & marg_array(:,2)==sel(2) & marg_array(:,3)==sel(3)};
        G = (diag(lambda)*(A{ind(2)}'*A{ind(2)}).*(A{ind(1)}'*A{ind(1)})*diag(lambda));
        [~,mode] = find(sel==n);
        V = mtkrprod(X,A(sel),mode)*diag(lambda);
        grad=-2*V+2*A{n}*G;
        grad_array{n}(:,j)=grad(:);
    end
end
grad_array{N+1}=zeros(F,J(N+1));
for j=1:J(N+1)
    sel=marg_array(j,:);
    X=Y{marg_array(:,1)==sel(1) & marg_array(:,2)==sel(2) & marg_array(:,3)==sel(3)};
    V = khatrirao(A(sel),'r');
    G = ((A{sel(3)}'*A{sel(3)}).*(A{sel(2)}'*A{sel(2)}).*(A{sel(1)}'*A{sel(1)}));
    M = X(:);
    grad=-2*(V')*M+2*G*lambda;
    grad_array{N+1}(:,j)=grad(:);
end
        
end