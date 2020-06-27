
function long_x = get_triplets(X)
%get_triplets
%prepare the long_x vector, which consists of all vectorized marginal
%triplets PMFs
%
%By Arie Yeredor, Nov. 2018, revised May 2020
long_x=[];
for i=1:length(X)
    long_x=[long_x;X{i}(:)];
end
end
