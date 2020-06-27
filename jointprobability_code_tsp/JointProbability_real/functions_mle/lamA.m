function X=lamA(lam,A)
%calculate the full tensor from its factors
%lam is the loading vector of size F
%AA{N} is a cell array holding the factor matrices, AA{n} is In by F
%By Arie Yeredor, Nov. 2018, revised May 2020

%extract dimensions
F=length(lam);
N=length(A);
Ivec=zeros(1,N);
for n=1:N;
    Ivec(n)=size(A{n},1);
end
X=zeros(Ivec);
%calculate
for f=1:F
    xv=lam(f);
    for n=N:-1:1
        xv=kron(xv,A{n}(:,f));
    end
    X=X+reshape(xv,Ivec);
end