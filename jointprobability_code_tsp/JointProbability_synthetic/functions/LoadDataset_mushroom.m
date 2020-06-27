function [Data,I] = LoadDataset_mushroom()

fid = fopen('UCI_MUSHROOM.data');
A = textscan(fid,'%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c','delimiter',',');

c = size(A,2);
n_samples = size(A{1},1);
Data = zeros(n_samples,c);

% Poisonous 0/ Edible 1
Data(A{1}=='p',1)=0;
Data(A{1}=='e',1)=1;

% cap-shape: bell=b,conical=c,convex=x,flat=f, knobbed=k,sunken=s 
Data(A{2}=='b',2)=0;
Data(A{2}=='c',2)=1;
Data(A{2}=='x',2)=2;
Data(A{2}=='f',2)=3;
Data(A{2}=='k',2)=4;
Data(A{2}=='s',2)=5;

% cap-surface: fibrous=f,grooves=g,scaly=y,smooth=s 
Data(A{3}=='f',3)=0;
Data(A{3}=='g',3)=1;
Data(A{3}=='y',3)=2;
Data(A{3}=='s',3)=3;

% cap-color: brown=n,buff=b,cinnamon=c,gray=g,green=r, pink=p,purple=u,red=e,white=w,yellow=y 
Data(A{4}=='n',4)=0;
Data(A{4}=='b',4)=1;
Data(A{4}=='c',4)=2;
Data(A{4}=='g',4)=3;
Data(A{4}=='r',4)=4;
Data(A{4}=='p',4)=5;
Data(A{4}=='u',4)=6;
Data(A{4}=='e',4)=7;
Data(A{4}=='w',4)=8;
Data(A{4}=='y',4)=9;

% bruises: bruises=t,no=f 
Data(A{5}=='t',5)=0;
Data(A{5}=='f',5)=1;

% odor: almond=a,anise=l,creosote=c,fishy=y,foul=f, musty=m,none=n,pungent=p,spicy=s 
Data(A{6}=='a',6)=0;
Data(A{6}=='l',6)=1;
Data(A{6}=='c',6)=2;
Data(A{6}=='y',6)=3;
Data(A{6}=='f',6)=4;
Data(A{6}=='m',6)=5;
Data(A{6}=='n',6)=6;
Data(A{6}=='p',6)=7;
Data(A{6}=='s',6)=8;

% gill-attachment: attached=a,descending=d,free=f,notched=n 
Data(A{7}=='a',7)=0;
Data(A{7}=='f',7)=1;

% gill-spacing: close=c,crowded=w,distant=d 
Data(A{8}=='c',8)=0;
Data(A{8}=='w',8)=1;

% gill-size: broad=b,narrow=n 
Data(A{9}=='b',9)=0;
Data(A{9}=='n',9)=1;

% gill-color: blck=k,brn=n,bf=b,choco=h,gray=g,green=r,orange=o,pink=p,purple=u,red=e,white=w,yellow=y 
Data(A{10}=='k',10)=0;
Data(A{10}=='n',10)=1;
Data(A{10}=='b',10)=2;
Data(A{10}=='h',10)=3;
Data(A{10}=='g',10)=4;
Data(A{10}=='r',10)=5;
Data(A{10}=='o',10)=6;
Data(A{10}=='p',10)=7;
Data(A{10}=='u',10)=8;
Data(A{10}=='e',10)=9;
Data(A{10}=='w',10)=10;
Data(A{10}=='y',10)=11;

% stalk-shape: enlarging=e,tapering=t 
Data(A{11}=='e',11)=0;
Data(A{11}=='t',11)=1;

% stalk-root: bulbous=b,club=c,cup=u,equal=e,rhizomorphs=z,rooted=r,missing=?  
Data(A{12}=='b',12)=0;
Data(A{12}=='c',12)=1;
Data(A{12}=='e',12)=2;
Data(A{12}=='r',12)=3;
Data(A{12}=='?',12)=-1;

% stalk-surface-above-ring: fibrous=f,scaly=y,silky=k,smooth=s 
Data(A{13}=='f',13)=0;
Data(A{13}=='y',13)=1;
Data(A{13}=='k',13)=2;
Data(A{13}=='s',13)=3;

% stalk-surface-below-ring: fibrous=f,scaly=y,silky=k,smooth=s 
Data(A{14}=='f',14)=0;
Data(A{14}=='y',14)=1;
Data(A{14}=='k',14)=2;
Data(A{14}=='s',14)=3;

% stalk-color-above-ring: brown=n,buff=b,cinnamon=c,gray=g,orange=o, pink=p,red=e,white=w,yellow=y 
Data(A{15}=='n',15)=0;
Data(A{15}=='b',15)=1;
Data(A{15}=='c',15)=2;
Data(A{15}=='g',15)=3;
Data(A{15}=='o',15)=4;
Data(A{15}=='p',15)=5;
Data(A{15}=='e',15)=6;
Data(A{15}=='w',15)=7;
Data(A{15}=='y',15)=8;

% stalk-color-below-ring: brown=n,buff=b,cinnamon=c,gray=g,orange=o, pink=p,red=e,white=w,yellow=y 
Data(A{16}=='n',16)=0;
Data(A{16}=='b',16)=1;
Data(A{16}=='c',16)=2;
Data(A{16}=='g',16)=3;
Data(A{16}=='o',16)=4;
Data(A{16}=='p',16)=5;
Data(A{16}=='e',16)=6;
Data(A{16}=='w',16)=7;
Data(A{16}=='y',16)=8;

% veil-type: partial=p,universal=u 
% Data(A{17}=='p',17)=0;
% Data(A{17}=='u',17)=1;

% veil-color: brown=n,orange=o,white=w,yellow=y 
Data(A{18}=='n',18)=0;
Data(A{18}=='o',18)=1;
Data(A{18}=='w',18)=2;
Data(A{18}=='y',18)=3;

% ring-number: none=n,one=o,two=t 
Data(A{19}=='n',19)=0;
Data(A{19}=='o',19)=1;
Data(A{19}=='t',19)=2;

% ring-type: cobwebby=c,evanescent=e,flaring=f,large=l,none=n,pendant=p,sheathing=s,zone=z 
Data(A{20}=='e',20)=0;
Data(A{20}=='f',20)=1;
Data(A{20}=='l',20)=2;
Data(A{20}=='n',20)=3;
Data(A{20}=='p',20)=4;
Data(A{20}=='s',20)=5;

% spore-print-color: black=k,brown=n,buff=b,chocolate=h,green=r,orange=o,purple=u,white=w,yellow=y 
Data(A{21}=='k',21)=0;
Data(A{21}=='n',21)=1;
Data(A{21}=='b',21)=2;
Data(A{21}=='h',21)=3;
Data(A{21}=='r',21)=4;
Data(A{21}=='o',21)=5;
Data(A{21}=='u',21)=6;
Data(A{21}=='w',21)=7;
Data(A{21}=='y',21)=8;

% population: abundant=a,clustered=c,numerous=n, scattered=s,several=v,solitary=y
Data(A{22}=='a',22)=0;
Data(A{22}=='c',22)=1;
Data(A{22}=='n',22)=2;
Data(A{22}=='s',22)=3;
Data(A{22}=='v',22)=4;
Data(A{22}=='y',22)=5;

% habitat: grasses=g,leaves=l,meadows=m,paths=p,urban=u,waste=w,woods=d
Data(A{23}=='g',23)=0;
Data(A{23}=='l',23)=1;
Data(A{23}=='m',23)=2;
Data(A{23}=='p',23)=3;
Data(A{23}=='u',23)=4;
Data(A{23}=='w',23)=5;
Data(A{23}=='d',23)=6;

most_common = mode(Data);
Data(Data(:,12)==-1,12) = most_common(12);

Data= Data(:,[1:16 18:end]);

% number of different values each variable can take
I = zeros(1,c-1);
for i=1:c-1
    I(i) = length(unique(Data(:,i)));
end
% replace missing value with most common
Data = Data+1;
end
