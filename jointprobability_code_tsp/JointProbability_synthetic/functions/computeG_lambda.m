function G = computeG_lambda(marg,GG,F)
G = zeros(F,F);
for i=1:size(marg,1)
    G_ = ones(F,F);
    for j = marg{i}
        G_ = G_ .* GG{j};
    end
    G = G + G_;
end