function G = computeG(marg,rows,cols,lambda,GG)
F = length(lambda); 
G = zeros(F,F);

for i=1:length(rows)
    G_ = ones(F,F);
    for j =  marg{rows(i)}([1:cols(i)-1 cols(i)+1:end])
        G_ = G_ .* GG{j};
    end
    G = G + G_;
end
G = (lambda*lambda').*G;
