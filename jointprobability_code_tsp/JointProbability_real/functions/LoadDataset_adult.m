function [Data,I] = LoadDataset_adult()

fid = fopen('UCI_ADULT.data');
A = textscan(fid,'%f%s%f%s%f%s%s%s%s%s%f%f%f%s%s','delimiter',',');

% fid = fopen('UCI_ADULT2.data');
% B = textscan(fid,'%f%s%f%s%f%s%s%s%s%s%f%f%f%s%s','delimiter',',');
% 
% for i= 1:size(A,2)
%     A{i}=[A{i};B{i}];
% end
n_samples = length(A{1});
n_features = size(A,2);
Data = zeros(n_samples,n_features);

%%%%%%First attribute is age%%%%%%%%
min1 = min(A{1});
max1 = max(A{1});
step1 = 3;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{1}>=range(k) & A{1} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,1)=count;
    end
end
Data(A{1}>=range(length(range)),1)=count+1;


%%%%%%Second attribute is work class%%%%%
range = ["Private","Self-emp-not-inc","Self-emp-inc","Federal-gov","Local-gov","State-gov","Without-pay","Never-worked"];
for k= 1:length(range)
    Data(strcmp(A{2},range(k)),2)=k;
end

%%%% Third attribute is fnlwgt%%%%%%
min1 = min(A{3});
max1 = max(A{3});
step1 = 50000;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{3}>=range(k) & A{3} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,3)=count;
    end
end
Data(A{3}>=range(length(range)),3)=count+1;

%%%%% Fourth attribute is education%%%%
range = ["Bachelors", "Some-college", "11th", "HS-grad", "Prof-school", "Assoc-acdm", "Assoc-voc", "9th", "7th-8th", "12th", "Masters", "1st-4th", "10th", "Doctorate", "5th-6th", "Preschool"];
for k= 1:length(range)
    Data(strcmp(A{4},range(k)),4)=k;
end

%%%%%%Fifth attribute is education num%%%%%%%%
min1 = min(A{5});
max1 = max(A{5});
step1 = 1;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{5}>=range(k) & A{5} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,5)=count;
    end
end
Data(A{5}>=range(length(range)),5)=count+1;

%%%%% Sixth attribute is marital status%%%%
range = ["Married-civ-spouse", "Divorced", "Never-married", "Separated", "Widowed", "Married-spouse-absent", "Married-AF-spouse"];
for k= 1:length(range)
    Data(strcmp(A{6},range(k)),6)=k;
end

%%%%% Seventh attribute is ocuupation%%%%
range = ["Tech-support", "Craft-repair", "Other-service", "Sales", "Exec-managerial", "Prof-specialty", "Handlers-cleaners", "Machine-op-inspct", "Adm-clerical", "Farming-fishing", "Transport-moving", "Priv-house-serv", "Protective-serv", "Armed-Forces"];
for k= 1:length(range)
    Data(strcmp(A{7},range(k)),7)=k;
end

%%%%% Eighth attribute is relationship%%%%
range = ["Wife", "Own-child", "Husband", "Not-in-family", "Other-relative", "Unmarried"];
for k= 1:length(range)
    Data(strcmp(A{8},range(k)),8)=k;
end

%%%%% Ninth attribute is ethnicity%%%%
range =["White", "Asian-Pac-Islander", "Amer-Indian-Eskimo", "Other", "Black"];
for k= 1:length(range)
    Data(strcmp(A{9},range(k)),9)=k;
end

%%%%% Tenth attribute is sex%%%%
range = ["Female", "Male"];
for k= 1:length(range)
    Data(strcmp(A{10},range(k)),10)=k;
end

%%%%%%11th attribute is capital gain%%%%%%%%
min1 = min(A{11});
max1 = max(A{11});
step1 = 5000;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{11}>=range(k) & A{11} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,11)=count;
    end
end
Data(A{11}>=range(length(range)),11)=count+1;

%%%%%%12th attribute is capital loss%%%%%%%%
min1 = min(A{12});
max1 = max(A{12});
step1 = 1000;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{12}>=range(k) & A{12} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,12)=count;
    end
end
Data(A{12}>=range(length(range)),12)=count+1;

%%%%%%13th attribute is hours per week%%%%%%%%
min1 = min(A{13});
max1 = max(A{13});
step1 = 5;
range = min1:step1:max1;
count=0;
for k= 1:length(range)-1
    u=(A{13}>=range(k) & A{13} < range(k+1));
    if(sum(u)~=0)
        count=count+1;
        Data(u,13)=count;
    end
end
Data(A{13}>=range(length(range)),13)=count+1;

%%%%% 14th attribute is native country%%%%
range = ["United-States", "Cambodia", "England", "Puerto-Rico", "Canada", "Germany", "Outlying-US(Guam-USVI-etc)", "India", "Japan", "Greece", "South", "China", "Cuba", "Iran", "Honduras", "Philippines", "Italy", "Poland", "Jamaica", "Vietnam", "Mexico", "Portugal", "Ireland", "France", "Dominican-Republic", "Laos", "Ecuador", "Taiwan", "Haiti", "Columbia", "Hungary", "Guatemala", "Nicaragua", "Scotland", "Thailand", "Yugoslavia", "El-Salvador", "Trinadad&Tobago", "Peru", "Hong", "Holand-Netherlands"];
for k= 1:length(range)
    Data(strcmp(A{14},range(k)),14)=k;
end

%%%%% 14th attribute is income prediction%%%%
Data( strcmp(A{15},'>50K')    ,15)    = 1;
Data( strcmp(A{15},'>50K.')   ,15)    = 1;

Data( strcmp(A{15},'<=50K')   ,15)    = 2;
Data( strcmp(A{15},'<=50K.')  ,15)    = 2;

 %Remove rows having undefinied values
 for i =1:n_features
     Data(Data(:,i)==0,:)=[];
 end
 
 Data = Data(:,[n_features 2:n_features-1 1]);
    
I = zeros(1,n_features);
for i=1:n_features
    I(i) = length(unique(Data(:,i)));
end


end