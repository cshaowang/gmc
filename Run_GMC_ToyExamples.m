%% Toy examples on Two-Moon data set and Three-Ring data set
% Graph-based Multi-view Clustering (GMC)
% 
%%
clc;  close all; clear all;
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Load toy data set
datadir = 'Dataset/';
m = 2; % number of views
X = cell(1,m);
dataname = input('Input the name of data set: (TwoMoon or ThreeRing)\n','s');
while(1)
    if strcmp(dataname,'TwoMoon')
        dataf = [datadir, dataname];
        c = 2;
        load(dataf);
        break;
    elseif strcmp(dataname,'ThreeRing')
        dataf = [datadir, dataname];
        c = 3;
        load(dataf);
        flag = 1;
        break;
    else
        dataType = input('Please only input TwoMoon or ThreeRing\n','s');
    end;
    
end

%% Call GMC algorithm
num = size(X{1},1); % the number of samples
data = cell(1,m);
for i = 1:m
    data{i} = X{i}';
end
[predY, U, S0, S0_initial, F, evs] = GMC(data, c, 1, 0);
metric = CalcMeasures(y0(:,1), predY);
fprintf('Data set %s-> ACC:%.4f\tNMI:%.4f\tARI:%.4f\terror_cnt:%d\n',dataname,metric(1),metric(2),metric(3),metric(4));

markerSize = 20;
%% Original data
for v = 1:m
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', 20); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axix
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'The orighnal data';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% Original connected graph with probabilistic neighbors, line width denotes similarity
S1 = cell(1,m);
for v = 1:m
    S1{v} = S0_initial{v};
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = S1{v}(ii, jj);
            if weight > 0
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) %set y-axis
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Connected graph with probabilistic neighbors';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% the learned graph of each view by GMC, line width denotes similarity
S2 = cell(1,m);
for v = 1:m
    S2{v} = S0{v};
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = S2{v}(ii, jj);
            if weight > 0
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axis
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Connected graph with probabilistic neighbors';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% the learned unified graph by GMC, line width denotes similarity
U2 = U;
for v = 1:m
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = U2(ii, jj);
            if weight > 0
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axix
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Learnt connected graph';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;
