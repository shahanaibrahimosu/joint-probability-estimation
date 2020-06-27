function [Data,I] = LoadDataset_nursery()
fid = fopen('UCI_NURSERY.data');
A = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter',',');
dscr_features = [1:5 7:9];
c = length(dscr_features);

n_samples = size(A{1},1);
Data = zeros(n_samples,c);

Data( strcmp(A{1},'great_pret')          ,1) = 1;
Data( strcmp(A{1},'pretentious')          ,1) = 2;
Data( strcmp(A{1},'usual')          ,1) = 3;

Data( strcmp(A{2},'very_crit')          ,2) = 1;
Data( strcmp(A{2},'critical')          ,2) = 2;
Data( strcmp(A{2},'improper')          ,2) = 3;
Data( strcmp(A{2},'less_proper')          ,2) = 4;
Data( strcmp(A{2},'proper')          ,2) = 5;

Data( strcmp(A{3},'foster')        ,3) = 1;
Data( strcmp(A{3},'incomplete')    ,3) = 2;
Data( strcmp(A{3},'completed')     ,3) = 3;
Data( strcmp(A{3},'complete')      ,3) = 4;

Data( strcmp(A{4},'1')          ,4) = 1;
Data( strcmp(A{4},'2')          ,4) = 2;
Data( strcmp(A{4},'3')      ,4) = 3;
Data( strcmp(A{4},'more')      ,4) = 4;

Data( strcmp(A{5},'critical')      ,5) = 1;
Data( strcmp(A{5},'less_conv')      ,5) = 2;
Data( strcmp(A{5},'convenient')      ,5) = 3;

% Data( strcmp(A{6},'inconv')      ,6) = 1;
% Data( strcmp(A{6},'convenient')      ,6) = 2;

Data( strcmp(A{7},'problematic')    ,6) = 1;
Data( strcmp(A{7},'slightly_prob')      ,6) = 2;
Data( strcmp(A{7},'nonprob')     ,6) = 3;

Data( strcmp(A{8},'not_recom')    ,7) = 1;
Data( strcmp(A{8},'priority')      ,7) = 2;
Data( strcmp(A{8},'recommended')     ,7) = 3;

Data( strcmp(A{9},'not_recom')  ,8) = 1;
Data( strcmp(A{9},'recommend'),8) = 2;
Data( strcmp(A{9},'very_recom')     ,8) = 3;
Data( strcmp(A{9},'priority')     ,8) = 4;
Data( strcmp(A{9},'spec_prior')     ,8) = 4;

Data = Data(:,[c 1:c-1]);
I = zeros(1,c);

for i=1:c
    I(i) = length(unique(Data(:,i)));
end
end

