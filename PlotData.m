% Script to plot data %

%---Load Data---%
load('CorrDistributions2.mat','MaxCosSim');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');

%%
%---Figure 1a: Max Correlation---%
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.1], [0.05 0.5]);
for i=1:size(DataSet,1)
    subplot(3,3,i);
    cult = DataSet{i,1}.culture;
    ch = DataSet{i,1}.channel;
    ic = DataSet{i,1}.ic;
    title(['Culture ',num2str(cult),' Channel: ',num2str(ch)]);
    [orderedCmat,~,~]=DendrogramOrderMaxCosSim(MaxCosSim{i},ic);
    imagesc(orderedCmat);
    set(gca,'XTick',1:size(ic,2)/2:size(ic,2), 'XTickLabel',{'','Neuro',''});
    set(gca,'YTick',size(ic,2)+1:(size(MaxCosSim{i},1)-size(ic,2)+1)/2:size(MaxCosSim{i},1), 'YTickLabel',{'','Astro',''},'YTickLabelRotation',-90);
    set(gca,'PlotBoxAspectRatio',[1,1,1]);
end
%%
%---Figure 2: Projection of Correlation on Realspace---%
which=6;
[e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValues(which,500);
figure;
PlotResultsOnMEA([],[],[],mask,ch,a1,a_e2,aval);
figure;
loadcult(which);
PlotA2Arealspace(mask,a2a,1000);

%%
%---Figure 3: Joint PSTH---%
load('PSTHdata_byPeaks.mat');
%or
% load('PSTHdata_byStarts.mat');

figure;
i=6;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,4,1:2);
cult = DataSet{i,1}.culture;
ch = DataSet{i,1}.channel;
title(['Astrocyte Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
%---Order by Offsets---%
% order = SortPSTHpairs(atrigpsth{i});
%---No Ordering---%
order = 1:size(atrigpsth{i},1);
%----------------------%
temp = zscore(atrigpsth{i}(order,:)')';
temp(temp<0)=0;
temp = temp(2500:end,:); 
imagesc(atriglag{6},1:size(temp,1),temp);
xlim([-4900,4900]);
% set(gca,'XTick',[1,10:10:100],'XTickLabel',atriglag{i}(1):1000:atriglag{i}(end)+1000,'YTick',[]);
ylabel('Pairs');

subplot(2,4,3:4);
title(['Burst Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
temp = zscore(ntrigpsth{i}')';
temp(temp<0)=0;
imagesc(temp); 
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

%% 
%---Figure 4: Correlation vs. Distance---%
figure;
which=6;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
cr = matfile('CorrDist.mat');
a2as = cr.a2as(1,which); a2as=a2as{1};
a2ns = cr.a2ns(1,which);a2ns=a2ns{1};
n2ns = cr.n2ns(1,which);n2ns=n2ns{1};
a2ad = lg.a2ad(1,which); a2ad=a2ad{1};
a2nd = lg.a2nd(1,which);a2nd=a2nd{1};
n2nd = lg.n2nd(1,which);n2nd=n2nd{1};

subplot(3,1,1);
scatter(n2nD,n2nV);
title('Neuron to Neuron');
xlabel('Distance');
ylabel('Correlation');
subplot(3,1,2);
scatter(a2aD,a2aV);
title('Astrocyte to Astrocyte');
xlabel('Distance');
ylabel('Correlation');
subplot(3,1,3)
scatter(a2nD,a2nV);
title('Astrocyte to Neuron');
xlabel('Distance');
ylabel('Correlation');

%%
%---Figure 5: Delay vs. Distance---%
figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
which=6;
lg = matfile('LagDist2.mat');
a2aL = lg.a2aL(1,which); a2aL=a2aL{1};
a2nL = lg.a2nL(1,which);a2nL=a2nL{1};
n2nL = lg.n2nL(1,which);n2nL=n2nL{1};
a2ad = lg.a2ad(1,which); a2ad=a2ad{1};
a2nd = lg.a2nd(1,which);a2nd=a2nd{1};
n2nd = lg.n2nd(1,which);n2nd=n2nd{1};

subplot(3,1,1);
scatter(n2nd,abs(n2nL));
title('Neuron to Neuron');
xlabel('Distance');
ylabel('Delays at Max Corr');
subplot(3,1,2);
scatter(a2ad,abs(a2aL));
title('Astrocyte to Astrocyte');
xlabel('Distance');
ylabel('Delays at Max Corr');
subplot(3,1,3)
scatter(a2nd,abs(a2nL));
title('Astrocyte to Neuron');
xlabel('Distance');
ylabel('Delays at Max Corr');

%%
% Plot Distributions of A2A vs. A2N vs. N2N

clear all;
i=6;
loadcult(i);
a2n=a2n(:);
nn = tril(nan(size(n2n)),1)+triu(n2n,1);
n2n = nn(:);
n2n(isnan(n2n))=[];
aa = tril(nan(size(a2a)),1)+triu(a2a,1);
a2a = aa(:);
a2a(isnan(a2a))=[];
clear_all_but('a2a','a2n','n2n');

figure;
s1 = subplot(3,1,1);
s2 = subplot(3,1,2);
s3 = subplot(3,1,3);

dispHIST(a2n', OPTBINS(a2n',500));
set(gca,'XLim',[0,1]);
copyobj(allchild(gca),s1);


dispHIST(n2n', OPTBINS(n2n',500));
set(gca,'XLim',[0,1]);
copyobj(allchild(gca),s2);


dispHIST(a2a', OPTBINS(a2a',500));
set(gca,'XLim',[0,1]);
copyobj(allchild(gca),s3);

title(s1,'A2N');
title(s2,'N2N');
title(s3,'A2A');
