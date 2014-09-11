load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);

%%
%---Correlation vs. Distance---%
load('CorrDist.mat');
figure;
xlimit = 500;
ylimit = 1;
optM=[10,10];
for i=1:9
    subplot(3,3,i)
%     optM = a2aB{i};
%     a2a = [a2ad{i}(~isnan(a2as{i})),a2as{i}(~isnan(a2as{i}))]';
%     a2a = a2a(:,a2a(2,:)~=0);
    a2a = [a2ad{i},a2as{i}]';
    [~,bins]  = hist3([a2a(1,:)',a2a(2,:)'],'Edges',{linspace(min(a2a(1,:)),max(a2a(1,:)),optM(1)),linspace(0,max(a2a(2,:)),optM(2))});
    factor = Calc2DhistNormalization([],bins{1},[],~DataSet{i}.mask,'a2a');
    Plot2DHist(a2a(1,:)',a2a(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['A2A: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
xlimit = 4000;
optM=[10,35];
for i=1:9
    subplot(3,3,i)
%     optM = a2nB{i};
    a2n = [a2nd{i},a2ns{i}]';
    [~,bins]  = hist3([a2n(1,:)',a2n(2,:)'],'Edges',{linspace(min(a2n(1,:)),max(a2n(1,:)),optM(1)),linspace(0,max(a2n(2,:)),optM(2))});
    factor = Calc2DhistNormalization(DataSet{i}.channel,bins{1},MeaMap,size(a2n,2),'a2n');
    Plot2DHist(a2n(1,:)',a2n(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['A2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
optM=[100,10];
for i=1:9
    subplot(3,3,i)
%     optM = n2nB{1};
%     n2n = [n2nd{i}(~isnan(n2ns{i})),n2ns{i}(~isnan(n2ns{i}))]';
%     n2n = n2n(:,n2n(2,:)~=0);
    n2n = [n2nd{i},n2ns{i}]';
    [~,bins]  = hist3([n2n(1,:)',n2n(2,:)'],'Edges',{linspace(min(n2n(1,:)),max(n2n(1,:)),optM(1)),linspace(min(n2n(2,:)),max(n2n(2,:)),optM(2))});
    factor = Calc2DhistNormalization(DataSet{i}.ic(1,:),bins{1},MeaMap,[],'n2n');
    Plot2DHist(n2n(1,:)',n2n(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['N2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end
clear_all_but('DataSet');

%%
%---Lag versus Distance---%
clear_all_but('DataSet','subplot');
load('LagDist.mat');
load('MeaMapPlot.mat','MeaMap');

%--A2A--%
figure;
xlimit = 600;
ylimit = 500;
optM=[10,10];
for i=1:9
    subplot(3,3,i)
    %     optM = a2aB{i};
    a2a = [a2ad{i},abs(a2aL{i})]';
    [~,bins]  = hist3([a2a(1,:)',a2a(2,:)'],'Edges',{linspace(min(a2a(1,:)),max(a2a(1,:)),optM(1)),linspace(0,max(a2a(2,:)),optM(2))});
    factor = Calc2DhistNormalization([],bins{1},[],~DataSet{i}.mask,'a2a');
    Plot2DHist(a2a(1,:)',a2a(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['A2A: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

%--A2N--%
figure;
xlimit = 4000;
optM=[10,110];
for i=1:9
    subplot(3,3,i)
    %     optM = a2nB{i};
    a2n = [a2nd{i},abs(a2nL{i})]';
    [~,bins]  = hist3([a2n(1,:)',a2n(2,:)'],'Edges',{linspace(min(a2n(1,:)),max(a2n(1,:)),optM(1)),linspace(min(a2n(2,:)),max(a2n(2,:)),optM(2))});
    factor = Calc2DhistNormalization(DataSet{i}.channel,bins{1},MeaMap,size(a2n,2),'a2n');
    Plot2DHist(a2n(1,:)',a2n(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['A2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

%--N2N--%
figure;
optM=[100,100];
for i=1:9
    subplot(3,3,i)
    %     optM = n2nB{1};
    n2n = [n2nd{i},abs(n2nL{i})]';
    [~,bins]  = hist3([n2n(1,:)',n2n(2,:)'],'Edges',{linspace(min(n2n(1,:)),max(n2n(1,:)),optM(1)),linspace(min(n2n(2,:)),max(n2n(2,:)),optM(2))});
    factor = Calc2DhistNormalization(DataSet{i}.ic(1,:),bins{1},MeaMap,[],'n2n');
    Plot2DHist(n2n(1,:)',n2n(2,:)',optM(1),optM(2),xlimit,ylimit,factor,bins);
    title(['N2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

%%
%---Lag Distributions---%

% parfor i=1:9
% [aaHisL{i},aaBinL{i}]=CalcHist(a2aL{i});
% [anHisL{i},anBinL{i}]=CalcHist(a2nL{i});
% [nnHisL{i},nnBinL{i}]=CalcHist(n2nL{i});
% end
% save('LagHistos.mat');

figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
clear_all_but('DataSet');
load('LagHistos.mat');
c=reshape(1:27,3,9);
c=c';

for i=1:9
    subplot(3,9,c(i,1))
    offset =nanmin(aaHisL{i}(aaHisL{i}~=0));
    bar(aaBinL{i},log10(aaHisL{i}+offset),'histc');
    ylim([0,4]);
    yt=get(gca,'ytick');
    ytl = arrayfun(@(x) num2str(x),yt,'UniformOutput',false);
    ytl = cellfun(@(x) strcat('10^{',x,'}'),ytl,'UniformOutput',false);
    set(gca,'YTick',yt,'YTickLabel',ytl);
    xlim([-200,200]);
    xlabel('A2A');
    
    subplot(3,9,c(i,2))
    offset =nanmin(anHisL{i}(anHisL{i}~=0));
    bar(anBinL{i},log10(anHisL{i}+offset),'histc');
    ylim([0,4]);
    yt=get(gca,'ytick');
    ytl = arrayfun(@(x) num2str(x),yt,'UniformOutput',false);
    ytl = cellfun(@(x) strcat('10^{',x,'}'),ytl,'UniformOutput',false);
    set(gca,'YTick',yt,'YTickLabel',ytl);
    xlim([-200,200]);
    title(['Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
    xlabel('A2N');
    
    subplot(3,9,c(i,3))
    offset =nanmin(nnHisL{i}(nnHisL{i}~=0));
    bar(nnBinL{i},log10(nnHisL{i}+offset),'histc');
    ylim([0,4]);
    yt=get(gca,'ytick');
    ytl = arrayfun(@(x) num2str(x),yt,'UniformOutput',false);
    ytl = cellfun(@(x) strcat('10^{',x,'}'),ytl,'UniformOutput',false);
    set(gca,'YTick',yt,'YTickLabel',ytl);
    xlim([-200,200]);
    xlabel('N2N');
end

%%
%---Corr Distributions---%

% parfor i=1:9
% [aaHiss{i},aaBins{i}]=CalcHist(a2as{i});
% [anHiss{i},anBins{i}]=CalcHist(a2ns{i});
% [nnHiss{i},nnBins{i}]=CalcHist(n2ns{i});
% end
% save('CorrHistos.mat');

figure;
clear_all_but('DataSet');
load('CorrHistos.mat');
c=reshape(1:27,3,9);
c=c';

for i=1:9
    subplot(3,9,c(i,1))
    bar(aaBinS{i},aaHisS{i},'histc');
    %     xlim([0,1000]);
    xlabel('A2A');
    
    subplot(3,9,c(i,2))
    bar(anBinS{i},anHisS{i},'histc');
    %     xlim([0,1000]);
    title(['Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
    xlabel('A2N');
    
    subplot(3,9,c(i,3))
    bar(nnBinS{i},nnHisS{i},'histc');
    %     xlim([0,1000]);
    xlabel('N2N');
end
