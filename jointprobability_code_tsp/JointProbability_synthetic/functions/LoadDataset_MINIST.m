function [Data,I] = LoadDataset_MINIST()

load('Datasets/MINIST.mat'); % training data stored in arrays X, y
Data = [y X];
n_features = size(X,2);
I = zeros(1,c);
for i=1:c
    I(i) = length(unique(Data(:,i)));
end

end