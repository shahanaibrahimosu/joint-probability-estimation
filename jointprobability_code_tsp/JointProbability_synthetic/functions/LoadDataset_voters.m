function [Data,I] = LoadDataset_voters()
fid = fopen('UCI_VOTERS.data');
A = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');

n_samples = size(A{1},1);
c = size(A,2);

Data = zeros(n_samples,c);
for i=1:c
    Data(:,i) = A{i};
end

most_common = mode(Data);
for i = 1:c
    Data( Data(:,i)==2 ,i) = most_common(i);
end

% number of different values each variable can take
I = zeros(1,c);
for i=1:c
    I(i) = length(unique(Data(:,i)));
end
Data = Data + 1;
end