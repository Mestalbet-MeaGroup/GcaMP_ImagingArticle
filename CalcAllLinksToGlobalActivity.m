load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat');
load('CorrDistributions2.mat', 'MaxCosSim');
nnscores = [];
nascores = [];
anscores = [];
aascores = [];
nindex   = [];
aindex   = [];
cults    = [];
for i=1:9
    corrs = MaxCosSim{i};
    corrs = triu(corrs,1)+tril(nan(size(corrs)),1);
    numc = size(DataSet{i}.ic,2);
    n2n = corrs(1:numc,1:numc);
    a2a = corrs(numc+1:end,numc+1:end);
    a2n = corrs(1:numc,numc+1:end);
    
%     nnscores = [nnscores;nanmean(n2n,2)];
%     nascores = [nascores;nanmean(a2n,2)];
%     anscores = [anscores,nanmean(a2n,1)];
%     aascores = [aascores;nanmean(a2a,2)];
%     
    %     for j=1:size(n2n,1), numspks(j) = numel(DataSet{i}.t(DataSet{i}.ic(3,j):DataSet{i}.ic(4,j))); end
        nnscores = [nnscores;max(n2n,[],2)];
        nascores = [nascores;max(a2n,[],2)];
        [temp,ni] = max(a2n,[],1);
        anscores = [anscores,temp]; %pairwise correlation
        nindex = [nindex,ni];
        aindex = [aindex,1:numel(temp)];
        cults = [cults,ones(size(ni)).*i];
        aascores = [aascores;max(a2a,[],2)];   
end
clear n2n;
clear a2n;
clear a2a;

for i=1:9
[neuro{i},astro{i},n2a{i},a2n{i}] = CalcLinkToNetwork(DataSet,i);
end
% Neuro = N to GFR
% Astro = A to Mean Trace
% N2A = Neuron to mean trace
% A2N = Astrocyte to GFR
datatable=[];
for k=1:9
[Corr2Astro,BestAstro]=nanmax(MaxCosSim{k}(1:size(DataSet{k}.ic,2),size(DataSet{k}.ic,2)+1:end),[],2);
Astro2GFR = a2n{k}(BestAstro);
Channel = DataSet{k}.ic(1,:);
[Corr2GFR,ix] = sort(neuro{k},'descend');
temp = table(Channel(ix)',Corr2GFR',BestAstro(ix),Corr2Astro(ix),Astro2GFR','VariableNames',{'Channel','Corr2GFR','BestAstro','Corr2Astro','Astro2GFR'});
datatable=[datatable;temp];
end

%% Plot 1: Histograms
neuro1=cell2mat(neuro);
astro1=cell2mat(astro);
n2a1=cell2mat(n2a);
a2n1=cell2mat(a2n); %correlation to mean-field
%--Cult ID for each a2n value--%
% a2n1_cultID=[];
% sizes = cellfun(@(x) numel(x),a2n);
% for i=1:numel(sizes)
%     a2n1_cultID = [a2n1_cultID;ones(sizes(i),1).*i];
% end
%----%
figure;
bins=nsOPTBINS([neuro1,astro1,n2a1,a2n1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.005], [0.05 0.05], [0.05 0.05]);
subplot(4,1,1)
hold on;
[c,b]=hist(neuro1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'XTick',[]);
set(gca,'FontSize',18)
axis('tight');
xlim([0,1]);
title('Neurons to GFR');

subplot(4,1,2)
hold on;
[c,b]=hist(n2a1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'FontSize',18)
set(gca,'XTick',[]);
axis('tight');
xlim([0,1]);
title('Neurons to Average Astrocyte Trace');

subplot(4,1,3)
hold on;
[c,b]=hist(a2n1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'FontSize',18)
axis('tight');
xlim([0,1]);
set(gca,'XTick',[]);
title('Astrocytes to GFR');

subplot(4,1,4)
hold on;
[c,b]=hist(astro1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'FontSize',18)
axis('tight');
xlim([0,1]);
title('Astrocytes to Average Astrocyte Trace');

%% Plot 2: 2D Histograms
figure;
x = [nnscores;aascores]; y=[neuro1';a2n1']; 
z=[x,y];
clear x; 
clear y;
bins = nsOPTBINS(z(~isnan(z(:,1))&~isnan(z(:,2)),:)');
[~,bins]=hist3(z,bins);

subplot(1,3,1);
z = [nnscores,neuro1'];
hist3(z,'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
xlabel('Pairwise Neuron');
ylabel('Correlation to GFR');
set(gca,'FontSize',18);

subplot(1,3,2);
z = [aascores,a2n1'];
hist3(z,'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
set(gca,'YTick',[]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
xlabel('Pairwise Astrocyte');
set(gca,'FontSize',18);

subplot(1,3,3);
z = [nascores,neuro1'];
hist3(z,'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
set(gca,'YTick',[]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
xlabel('Pairwise Neuron to Astrocyte');
set(gca,'FontSize',18);

%% Plot 3:
z=[nascores;anscores'];
z=z(~isnan(z))';
bins = nsOPTBINS(z);

subplot(4,1,1);
[c,b]=hist(nascores,bins);
bar(b,c./trapz(b,c),'histc');
title('Max Pairwise Neuronal Correlations in A2N'); 
axis('tight');
xlim([0,0.6]);
set(gca,'XTick',[],'FontSize',18);

subplot(4,1,2)
hold on;
[c,b]=hist(n2a1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'FontSize',18)
set(gca,'XTick',[]);
axis('tight');
xlim([0,0.6]);
title('Neurons to Average Astrocyte Trace');

subplot(4,1,3);
[c,b]=hist(anscores,bins);
bar(b,c./trapz(b,c),'histc');
title('Max Pairwise Astrocyte Correlations in A2N');
set(gca,'FontSize',18);
axis('tight');
set(gca,'XTick',[]);
xlim([0,0.6]);

subplot(4,1,4)
hold on;
[c,b]=hist(a2n1,bins);
bar(b,c./trapz(b,c),'histc');
set(gca,'FontSize',18)
axis('tight');
xlim([0,0.6]);
title('Astrocytes to GFR');

%% Plot 4: Plot fit: Pairwise vs. MeanField
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.05 0.05]);
[xData, yData] = prepareCurveData( a2n1, anscores );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot residuals.
figure;
subplot(2,1,2);
h = plot( fitresult, xData, yData, 'stresiduals');
xdata = get(h,'XData'); % x-values
ydata = get(h,'YData'); %residual hight
hold on;
stem(xdata{1},ydata{1},...
    'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',5,...
    'Color',[0 0.447 0.741]);
stem(xdata{1}(find(ydata{1}>1.96)),ydata{1}(find(ydata{1}>1.96)),...
    'MarkerFaceColor',[237,47,89]./255,...
    'MarkerEdgeColor',[237,47,89]./255,...
    'MarkerSize',5,...
    'Color',[0 0.447 0.741]);
line(xdata{2},ydata{2});
line(xdata{2},ones(size(xdata{2}))*1.96,'LineStyle','--','Color','r');
ylabel('PairWise A2N [z-score]');
xlabel('Astrocyte to Neural Mean Field (average firing rate)');
axis tight;
grid on
set(gca,'TickDir','out');
set(gca,'FontSize',18);

winsize = cellfun(@(x) mean(diff(x.dfTime)),DataSet);
%Find the group of residuals these are; 
residuals = xdata{1}(find(ydata{1}>1.96));
BestAstros = find(ismember(a2n1,residuals));
BestNeuros(1,:) = nindex(BestAstros);
BestNeuros(2,:) = cults(BestAstros);
BestNeuros(3,:) = arrayfun(@(x,y) nanmean(DataSet{x}.FR(:,y))/winsize(x),BestNeuros(2,:),BestNeuros(1,:));
BestNeuros(4,:) = arrayfun(@(x,y) nanmean(MaxCosSim{x}(y,size(DataSet{x}.ic,2)+1:end)),BestNeuros(2,:),BestNeuros(1,:));
BestNeuros(5,:) = arrayfun(@(x,y) neuro{x}(y),BestNeuros(2,:),BestNeuros(1,:));
BestNeuros(6,:) = anscores(BestAstros);
BestNeuros(7,:) = aindex(BestAstros);
% HIBastros = intersect(find(a2n1>=mean(a2n1)+std(a2n1)),find(yData>=mean(yData)+2*std(yData)));
HIBastros = intersect(find(a2n1>=0.3),find(yData>=0.4)); %high in both

% HIBastros = find(ismember(a2n1,highinboth));
HIBn(1,:) = nindex(HIBastros);
HIBn(2,:) = cults(HIBastros);
HIBn(3,:) = arrayfun(@(x,y) nanmean(DataSet{x}.FR(:,y))/winsize(x),HIBn(2,:),HIBn(1,:));
HIBn(4,:) = arrayfun(@(x,y) nanmean(MaxCosSim{x}(y,size(DataSet{x}.ic,2)+1:end)),HIBn(2,:),HIBn(1,:));
HIBn(5,:) = arrayfun(@(x,y) neuro{x}(y),HIBn(2,:),HIBn(1,:));
HIBn(6,:) = anscores(HIBastros);
HIBn(7,:) = aindex(HIBastros);
otherFR     =  arrayfun(@(x) nanmean(nanmean(DataSet{x}.FR(:,setdiff(1:size(DataSet{x}.FR,2),unique(BestNeuros(1,BestNeuros(2,:)==x))))))/winsize(x),1:9);
otherCorr2A =  arrayfun(@(x) nanmean(nanmean(MaxCosSim{x}(setdiff(1:size(DataSet{x}.FR,2),unique(BestNeuros(1,BestNeuros(2,:)==x))),size(DataSet{x}.ic,2)+1:end))),1:9);
otherCorr2GFR = arrayfun(@(x,y) nanmean(neuro{x}(setdiff(1:size(DataSet{x}.FR,2),unique(BestNeuros(1,BestNeuros(2,:)==x))))),1:9);

subplot(2,1,1);
hold on
h = plot( fitresult, xData, yData);
plot(xData(BestAstros),yData(BestAstros),'.','Color',[237,47,89]./255); %color residuals above the confidence bound on the regression, red
set(h,'linewidth',5)
set(allchild(gca),'MarkerSize',10);
hold off;
ylabel('PairWise A2N [corr]');
axis tight;
grid on
set(gca,'XTickLabels',{},'TickDir','out');
xlabel('');
set(gca,'FontSize',18);


%---Figure of Neurons which are in high pairwise correlation with astro's but low correlation with the GFR---%
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
h = bar([1,2,3],[mean(BestNeuros(3,:)),mean(otherFR),mean(HIBn(3,:))],1);
errorbar([1,2,3],...
         [mean(BestNeuros(3,:)),mean(otherFR),mean(HIBn(3,:))],...
         [std(BestNeuros(3,:)),std(otherFR),std(HIBn(3,:))],'.k');
set(gca,'XTick',2,'XTickLabel',{'Mean FR'},'TickDir','out','YAxisLocation','left');
ylabel('<firing rate> [spikes*sec^-1]');
xlim([0,8]);
ylim([0,14]);
set(gca,'FontSize',18);

ha = axes('YAxisLocation','right','Color','none');
hold on;
h = bar(ha,[5,6,7],[mean(BestNeuros(5,:)),mean(otherCorr2GFR),mean(HIBn(5,:))],1);
errorbar([5,6,7],...
         [mean(BestNeuros(5,:)),mean(otherCorr2GFR),mean(HIBn(5,:))],...
         [std(BestNeuros(5,:)),std(otherCorr2GFR),std(HIBn(5,:))],'.k');
xlim([0,8]);
set(gca,'XTick',6,'XTickLabel',{'cross-correlation to GFR'},'TickDir','out');
ylabel(ha,'normalized cross-correlation'); 
set(gca,'FontSize',18);
% export_fig('Fig_ResidualsCorrFR.eps','-r600');

%% Statistics
g1 = [ones(size(BestNeuros(3,:))),ones(size(otherFR))*2,ones(size(HIBn(3,:)))*3]; %high res, all others, high both
g2 = [BestNeuros(2,:),1:9,HIBn(2,:)]; %which culture
frs = [BestNeuros(3,:),otherFR,HIBn(3,:)];
[~,~,stats] = anovan(frs,{g1,g2});
[frStat,~,~] = multcompare(stats,'Dimension',1);
corrs = [BestNeuros(5,:),otherCorr2GFR,HIBn(5,:)];
[~,~,stats] = anovan(corrs,{g1,g2});
[corrStat,~,~] = multcompare(stats,'Dimension',1);
