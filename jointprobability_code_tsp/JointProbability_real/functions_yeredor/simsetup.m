%simsetup
%simulation setup parameters
%
%By Arie Yeredor, Nov. 2018, revised May 2020


Ivec=[2 3 4 3 2];   %alphabet size for each variable
N=length(Ivec);     %number of variables
F=2;                %number of factors

lam=[0.7;0.3];

Aall=cell(1,N);
Aall{1}=[0.1 0.9;0.4 0.6]';
Aall{2}=[0.5 0.1 0.4;0.5 0.4 0.1]';
Aall{3}=[0.2 0.5 0.1 0.2;0.4 0.1 0.1 0.4]';
Aall{4}=[0.1 0.1 0.8;0.2 0.4 0.4]';
Aall{5}=[0.5 0.5;0.2 0.8]';

T=50000;   %observations length
