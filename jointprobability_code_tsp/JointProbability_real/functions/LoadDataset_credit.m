function [Data,I] = LoadDataset_credit()
fid = fopen('./Datasets/UCI_CREDIT.data');
A = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',' ');
dscr_features = [1 4 5 6 8 9 11 12 15];
c = length(dscr_features);

n_samples = size(A{1},1);

Data = zeros(n_samples,c);
for i=1:length(dscr_features)
    Data(:,i) = A{dscr_features(i)};
end

Data(:,2)  = Data(:,2) - 1; 
Data(:,3)  = Data(:,3) - 1; 
Data(:,4)  = Data(:,4) - 1; 
Data(:,8)  = Data(:,8) - 1; 

% move label to 1st position
Data = Data(:,[9 1:8]);
Data(Data(:,5)==8,5) = 5;

Data = Data + 1;
% number of different values each variable can take
I = zeros(1,c);
for i=1:c
    I(i) = length(unique(Data(:,i)));
end

end

