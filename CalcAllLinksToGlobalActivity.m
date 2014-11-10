load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat');
load('CorrDistributions2.mat', 'MaxCosSim');
nnscores = [];
nascores = [];
anscores = [];
aascores = [];
nindex = [];
cults = [];
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
        anscores = [anscores,temp];
        nindex = [nindex,ni];
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
a2n1=cell2mat(a2n);
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
subplot(2,1,1);
h = plot( fitresult, xData, yData);
set(h,'linewidth',5)
set(allchild(gca),'MarkerSize',10);
ylabel('PairWise A2N [corr]');
axis tight;
grid on
set(gca,'XTickLabels',{},'TickDir','out');
xlabel('');
set(gca,'FontSize',18);

subplot(2,1,2);
h = plot( fitresult, xData, yData, 'stresiduals');
xdata = get(h,'XData');
ydata = get(h,'YData');
hold on;
stem(xdata{1},ydata{1},...
    'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',5,...
    'Color',[0 0.447 0.741]);
line(xdata{2},ydata{2});
ylabel('PairWise A2N [z-score]');
xlabel('Astrocyte to Neural Mean Field (average firing rate)');
axis tight;
grid on
set(gca,'TickDir','out');
set(gca,'FontSize',18);

% Now brush all the positive residuals and copy them into a N*2 sized array
% named residual.

%Find the group of residuals these are; 

BestAstros = find(ismember(a2n1,residuals(find(residuals(:,2)>=2),1)));
BestNeuros(1,:) = nindex(BestAstros);
BestNeuros(2,:) = cults(BestAstros);
HIBastros = find(ismember(a2n1,highinboth(:,1)));
HIBn(1,:) = nindex(HIBastros);
HIBn(2,:) = cults(HIBastros);
HIBn(3,:) =  arrayfun(@(x,y) nanmean(DataSet{x}.FR(:,y)),HIBn(2,:),HIBn(1,:));

LGastros = find(ismember(a2n1,higha2nlowa2gfr(:,1)));
LGn(1,:) = nindex(LGastros);
LGn(2,:) = cults(LGastros);
LGn(3,:) =  arrayfun(@(x,y) nanmean(DataSet{x}.FR(:,y)),LGn(2,:),LGn(1,:));



BestNeuros(3,:) =  arrayfun(@(x,y) nanmean(DataSet{x}.FR(:,y)),BestNeuros(2,:),BestNeuros(1,:));
BestNeuros(4,:) =  arrayfun(@(x,y) nanmean(MaxCosSim{x}(y,size(DataSet{x}.ic,2)+1:end)),BestNeuros(2,:),BestNeuros(1,:));
BestNeuros(5,:) = arrayfun(@(x,y) neuro{x}(y),BestNeuros(2,:),BestNeuros(1,:));

% otherFR =  arrayfun(@(x,y) nanmean(nanmax(DataSet{x}.FR(:,setdiff(1:size(DataSet{x}.FR,2),y)))),BestNeuros(2,:),BestNeuros(1,:));
% otherCorr2A=  arrayfun(@(x,y) nanmean(nanmax(MaxCosSim{x}(setdiff(1:size(DataSet{x}.FR,2),y),size(DataSet{x}.ic,2)+1:end))),BestNeuros(2,:),BestNeuros(1,:));
otherFR =  arrayfun(@(x,y) nanmean(nanmean(DataSet{x}.FR(:,setdiff(1:size(DataSet{x}.FR,2),y)))),BestNeuros(2,:),BestNeuros(1,:));
otherCorr2A=  arrayfun(@(x,y) nanmean(nanmean(MaxCosSim{x}(setdiff(1:size(DataSet{x}.FR,2),y),size(DataSet{x}.ic,2)+1:end))),BestNeuros(2,:),BestNeuros(1,:));
otherCorr2GFR = arrayfun(@(x,y) nanmean(neuro{x}(setdiff(1:size(DataSet{x}.FR,2),y))),BestNeuros(2,:),BestNeuros(1,:));

figure;
hold on;
h = bar([mean(BestNeuros(3,:)),mean(otherFR);mean(BestNeuros(4,:)),mean(otherCorr2A);mean(BestNeuros(5,:)),mean(otherCorr2GFR)]);
w = get(h(1),'BarWidth');
bx = get(h(1),'XData');
errorbar([bx(1)-w/(numel(bx)*2),bx(1)+w/(numel(bx)*2);bx(2)-w/(numel(bx)*2),bx(2)+w/(numel(bx)*2);bx(3)-w/(numel(bx)*2),bx(3)+w/(numel(bx)*2)],...
         [mean(BestNeuros(3,:)),mean(otherFR);mean(BestNeuros(4,:)),mean(otherCorr2A);mean(BestNeuros(5,:)),mean(otherCorr2GFR)],...
         [std(BestNeuros(3,:)),std(otherFR);std(BestNeuros(4,:)),std(otherCorr2A);std(BestNeuros(5,:)),std(otherCorr2GFR)],'.')
set(gca,'XTickLabel',{'','Mean FR','','Max Corr to Astro','','Corr to GFR'});
legend({'High Residuals','All Others'});
set(gca,'FontSize',18);
axis tight;
xlim([0.5,3.5]);
