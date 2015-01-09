% Creates two panels of correlation figure:
% 1\Correlation of neurons and astrocytes within view-field  plotted on real space
% 2\ Histogram of correlations versus distance + distance vs. lag


%% Real space plot
howmany=500;
[e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValuesWithVF(6,howmany);
e1=e1(howmany+1:end);
e2=e2(howmany+1:end);
eval = eval(howmany+1:end);
ind = find(eval<min(aval));
e1(ind)=[];
e2(ind)=[];
eval(ind)=[];
PlotResultsOnMEAwithVF(e1,e2,eval,mask,ch,a1,a_e2,aval);



%% Plot Histograms
load('CorrDist2.mat');
d = cell2mat(a2nd');
s = abs(cell2mat(a2ns'));
if size(d,1)>size(d,2)
    d=d';
    s=s';
end
ds = [d;s]';
ds(isnan(ds(:,2)),:)=[];
clear_all_but('ds');
% optM = nsOPTBINS(ds');
hist3(ds,[15,15]);
az=0;%35
el=90;%32
view([az,el]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
CData = get(get(gca,'child'),'CData');
CData(CData<=0) = NaN;
CData(CData==Inf) = NaN;
% Find min/max
minC = min(CData(:));
maxC = max(CData(:));
% Convert data to log space
CData = log10(CData);
% Normalize the log scaled CData to have the same range as the ZData so
% that the colorbar scale will be correct
CData = minC + (maxC-minC)*(CData-log10(minC))/(log10(maxC/minC));
% Now set the CData of the surface to the normalized CData
set(get(gca,'child'),'CData',CData)
set(get(gca,'child'),'CDataMapping','scaled')
axis tight;
grid off;
set(gca,'PlotBoxAspectRatio',[1,1,1],'TickDir','out','FontSize',9);
ylabel('astrocyte-neuron pairwise correlation');
xlabel('distance [um]');
hCbar = colorbar(gca);
myscale = [nanmin(CData(:)),10^2,10^2.5,10^3,10^3.5,nanmax(CData(:))];
caxis([myscale(1),myscale(end)])
set(hCbar,'YTick',myscale);
set(hCbar,'YTickLabel',log10(myscale),'TickDir','out','FontSize',9);
ylabel(hCbar,'incidence [counts]');

%--------------------------------------Now for Lags-----------------------%
load('E:\LagDist3.mat');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')
load('MeaMapPlot.mat','MeaMap');
d = cell2mat(a2nd');
l = abs(cell2mat(a2nL'));
if size(d,1)>size(d,2)
    d=d';
    l=l';
end
dl = [d;l]';
dl(isnan(dl(:,2)),:)=[];
bin = {linspace(min(dl(:,1)),max(dl(:,1)),15),linspace(min(dl(:,2)),max(dl(:,2)),16)};
normhist=zeros(numel(bin{1}),numel(bin{2}));
DistMat = CreateElectrodeDistanceTable();
for i=1:size(DataSet,1)
    [hist,bins]  = hist3([a2nd{i}(~isnan(a2nd{i})),a2nL{i}(~isnan(a2nd{i}))],'Edges',bin);
    nnB = size(bins{1},2);
    factor=zeros(1,nnB(1));
    for k=1:nnB(1)-1
        factor(k) = NormalizeDistanceCounts(DistMat,MeaMap,DataSet{i}.channel,bins{1}(k),bins{1}(k+1));
    end
%     factor = factor./trapz(bins{1},factor);
    factor(factor==0)=inf;
    temphist = hist./repmat(factor,size(normhist,2),1)';
    normhist=normhist+temphist/trapz(bins{2},trapz(bins{1},temphist));
end
normhist = normhist./9;

xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
s = surface(xo,yo,normhist);
az=0;%35
el=90;%32
view([az,el]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
grid off;
hCbar = colorbar(gca);
set(gca,'PlotBoxAspectRatio',[1,1,1],'TickDir','out','FontSize',9);
ylabel('astrocyte-neuron lag at maximum correlation [ms]');
xlabel('distance [um]');
ylabel(hCbar,'incidence [counts]');
set(hCbar,'TickDirection','out');
