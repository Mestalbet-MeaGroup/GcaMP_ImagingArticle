% Script to plot data %

%---Load Data---%
CS  = matfile('F:\CosSimTemp.mat');
% db = matfile('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
load('CorrDistributions.mat','MaxCosSim');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
load('PSTHdata.mat')

%---Figure 1: Max Correlation---%
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.1], [0.05 0.5]);
for i=1:size(DataSet,1)
    subplot(3,3,i);
    cult = DataSet{i,1}.culture;
    ch = DataSet{i,1}.channel;
    ic = DataSet{i,1}.ic;
    title(['Culture ',num2str(cult),' Channel: ',num2str(ch)]);
    imagesc(MaxCosSim{i});
    set(gca,'XTick',1:size(ic,2)/2:size(ic,2), 'XTickLabel',{'','Neuro',''});
    set(gca,'YTick',size(ic,2)+1:(size(MaxCosSim{i},1)-size(ic,2)+1)/2:size(MaxCosSim{i},1), 'YTickLabel',{'','Astro',''},'YTickLabelRotation',-90);
    set(gca,'PlotBoxAspectRatio',[1,1,1]);
end

%---Figure 2: Projection of Correlation on Realspace---%
[e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValues(6,500);
figure;
PlotResultsOnMEA([],[],[],mask,ch,a1,a_e2,aval);
figure;
loadcult(6);
PlotA2Arealspace(mask,a2a,1000);


%---Figure 3: Joint PSTH---%
figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,4,1:2);
i=6;
cult = DataSet{i,1}.culture;
ch = DataSet{i,1}.channel;
ic = DataSet{i,1}.ic;
title(['Astrocyte Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
imagesc(atrigpsth{i});
xlim([1,100]);
set(gca,'XTick',[1,10:10:100],'XTickLabel',ntriglag{i}(1):1000:ntriglag{i}(end)+1000,'YTick',[]);
ylabel('Pairs');

subplot(2,4,3:4);
title(['Burst Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
i=6;
imagesc(ntrigpsth{i});
xlim([1,100]);
set(gca,'XTick',[1,10:10:100],'XTickLabel',ntriglag{i}(1):1000:ntriglag{i}(end)+1000,'YTick',[]);
ylabel('Pairs');

subplot(2,4,5:6)
hold on;
title('Astro-Triggered PSTH, sum over All Pairs');
for i=1:size(DataSet,1)
    plot(atriglag{i},nansum(zscore(atrigpsth{i},0,2),1),'.-');
end
subplot(2,4,7:8)
hold on;
title('Burst-Triggered PSTH, sum over All Pairs');
for i=1:size(DataSet,1)
    plot(ntriglag{i},nansum(zscore(ntrigpsth{i},0,2),1),'.-');
end

%---Figure 4: Correlation vs. Distance---%
% loadcult(6);
% [n2nD,n2nV] = CalcValueDist(n2n,ic,[],[],'n2n');
% [a2aD,a2aV] = CalcValueDist(a2a,[],mask,[],'a2a');
% [a2nD,a2nV] = CalcValueDist(a2n,ic,mask,ch,'a2n');
load('DistData.mat'); %continue here
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(1,3,1);
scatter(n2nD,n2nV);
title('Neuron to Neuron');
xlabel('Distance');
ylabel('Correlation');
subplot(1,3,2);
scatter(a2aD,a2aV);
title('Astrocyte to Astrocyte');
xlabel('Distance');
ylabel('Correlation');
subplot(1,3,3)
scatter(a2nD,a2nV);
title('Astrocyte to Neuron');
xlabel('Distance');
ylabel('Correlation');

%---Figure 5: Delay vs. Distance---%
figure;
[n2nL,a2aL,a2nL,lagmat,ic]=CalcLagsAtMaxCorr(6);
[n2nD,n2nV] = CalcValueDist(n2nL,ic,[],[],'n2n');
[a2aD,a2aV] = CalcValueDist(a2aL,[],mask,[],'a2a');
[a2nD,a2nV] = CalcValueDist(a2nL,ic,mask,ch,'a2n');
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(1,3,1);
scatter(n2nD,n2nV);
title('Neuron to Neuron');
xlabel('Distance');
ylabel('Lag at Max Corr');
subplot(1,3,2);
scatter(a2aD,a2aV);
title('Astrocyte to Astrocyte');
xlabel('Distance');
ylabel('Lag at Max Corr');
subplot(1,3,3)
scatter(a2nD,a2nV);
title('Astrocyte to Neuron');
xlabel('Distance');
ylabel('Lag at Max Corr');


%---Figure 6: Heatmaps---%
