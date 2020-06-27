%get_triplets
%prepare the long_x vector, which consists of all vectorized marginal
%triplets PMFs
%
%By Arie Yeredor, Nov. 2018, revised May 2020

long_x=[];
for j=1:N-2
    for k=j+1:N-1
        for l=k+1:N
            allix=1:N;
            allix([j k l])=[];
            Xjkl=Xemp;
            for nn=1:N-3;
                Xjkl=sum(Xjkl,allix(nn));
            end
            Xjkl=squeeze(Xjkl);
            long_x=[long_x;Xjkl(:)];
        end
    end
end
