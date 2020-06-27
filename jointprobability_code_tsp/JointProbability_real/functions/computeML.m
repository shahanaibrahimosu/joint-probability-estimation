function A_ml= computeML(dataset,F_true,I)
N=length(I);
A_ml     = cell(N,1);
ns=length(dataset);
count_F=zeros(1,F_true);
for n = 1 : N
    for i= 1: I(n)
        for f=1:F_true
            count_F(f)=sum(dataset(:,N+1)==f);
            A_ml{n}(i,f)=sum((dataset(:,n)==i) & (dataset(:,N+1)==f))/count_F(f);
        end
    end
    A_ml{n} = A_ml{n}*diag(1./sum(A_ml{n},1));
end
