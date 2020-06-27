function [X_comb,list_f] = construct_Z_Matrix(Y,I,marg,M)
marg_array=cell2mat(marg);
N=length(I);
K=length(Y{1});
X_comb = cell(N,1);
list_f = cell(N,1);
for i=1:N
    X=[];
    list=[];
    for j=[1:i-1 i+1:N]
        [ind,~]=sort([i,j]);
        X_sel = Y{marg_array(:,1)==ind(1) & marg_array(:,2)==ind(2)};
        if(rank(X_sel)== K)
            if(ind(1)~=i)
               X_sel=X_sel';
            end
            X = [X X_sel];
            list =[list j];
        end
    end
    X_comb{i}=X;
    list_f{i} = list;
end
end