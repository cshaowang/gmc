%% Experiments on real-world data sets
% Graph-based Multi-view Clustering (GMC)
%
%%
clc;  close all; clear all;
currentFolder = pwd;
addpath(genpath(currentFolder));
resultdir = 'Results/';
if(~exist('Results','file'))
    mkdir('Results');
    addpath(genpath('Results/'));
end
dataname = {'100leaves','3sources','BBC','BBCSport','HW','HW2sources','NGs','WebKB'};
runtimes = 1; % run-times on each dataset, default: 1
numdata = length(dataname);

for cdata = 1:numdata
%% read dataset
idata = cdata;
datadir = 'Dataset/';
dataf = [datadir, cell2mat(dataname(idata))];
load(dataf);

X = data;
y0 = truelabel{1};
c = length(unique(truelabel{1}));
%% iteration ...
for rtimes = 1:runtimes
    [y, U, S0, S0_initial, F, evs] = GMC(X, c); % c: the # of clusters
metric = CalcMeasures(y0, y);
ACC(rtimes) = metric(1);
NMI(rtimes) = metric(2);
ARI(rtimes) = metric(3);
error_cnt(rtimes) = metric(4);
disp(char(dataname(idata)));
fprintf('=====In iteration %d=====\nACC:%.4f\tNMI:%.4f\tARI:%.4f\terror_cnt:%d\n',rtimes,metric(1),metric(2),metric(3),metric(4));
end;
    Result(1,:) = ACC;
    Result(2,:) = NMI;
    Result(3,:) = ARI;
    Result(4,1) = mean(ACC);
    Result(4,2) = mean(NMI);
    Result(4,3) = mean(ARI);
    Result(5,1) = std(ACC);
    Result(5,2) = std(NMI);
    Result(5,3) = std(ARI);
save([resultdir,char(dataname(idata)),'_result.mat'],'Result','U','y0','y');
clear ACC NMI ARI metric Result U y0 y;
end;
