function [Data,I] = LoadDataset_heart()
fid = fopen('UCI_HEART.data');
A = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
dscr_features = [1 2 3 6 7 9 11 12 13 14];
c = length(dscr_features);

n_samples = size(A{1},1);
Data = zeros(n_samples,c);

for i=1:length(dscr_features)
    Data(:,i) = A{dscr_features(i)};
end

Data(:,3)   = Data(:,3)  - 1; 
Data(:,7)   = Data(:,7)  - 1; 

Data(A{13}==3,9)   = 0; 
Data(A{13}==6,9)   = 1; 
Data(A{13}==7,9)   = 2; 

Data(Data(:,10)>=1,10) = 1;

% move label to 1st position
Data = Data(:,[10 1:9]);

most_common = mode(Data);
for i = 1:c
    Data(Data(:,i)==-1,i) = most_common(i);
end

Data = Data(:,[1 3:end]);
Data = Data + 1;

c = 9;
% number of different values each variable can take
I = zeros(1,c);
for i=1:c
    I(i) = length(unique(Data(:,i)));
end
end

