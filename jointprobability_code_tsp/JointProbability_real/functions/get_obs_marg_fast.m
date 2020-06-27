function [Y,sum_y,len_list] = get_obs_marg_fast(Dataset,marg,I)
Y  = cell(size(marg,1),1);
sum_y  = zeros(size(marg,1),1);
len_list=cell(size(marg,1),1);
for i=1:size(marg,1)
    Y{i} = zeros(I(marg{i}));
    tmp = Dataset(all(Dataset(:,marg{i}) ~=0 ,2) , marg{i});
    [u,~,ic]=unique(tmp,'rows');
    C =[u histc(ic,1:size(u,1))];
    for c = 1:size(C,1)
        subscript = num2cell(C(c,1:end-1));
        ind = sub2ind(I(marg{i}),subscript{:});
        Y{i}(ind) = C(c,end);
    end
    len_list{i}=C(:,end);
    sum_y(i) = sum(Y{i}(:));
    Y{i} = Y{i} / sum(Y{i}(:));
    Y{i}(isnan(Y{i}))=0;
end
end