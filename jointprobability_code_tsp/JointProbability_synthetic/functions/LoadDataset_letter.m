function [Data,I] = LoadDataset_letter()
fid = fopen('letter-recognition.data');
A = textscan(fid,'%c%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
dscr_features = 1:17;
c = length(dscr_features);

n_samples = size(A{1},1);
Data = zeros(n_samples,c);

Data(:,1) = double(A{1})-64;
for i=2:length(dscr_features)
    Data(:,i) = A{dscr_features(i)};
end


% number of different values each variable can take
I = zeros(1,c);
for i=1:c
    I(i) = length(unique(Data(:,i)));
end
end

