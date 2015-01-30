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
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat')
PlotResultsOnMEAwithVF(e1,e2,eval,DataSet{6}.ic,mask,ch,a1,a_e2,aval);



%% Plot Histograms
load('CorrDist2.mat');
d = cell2mat(a2nd');
s = abs(cell2mat(a2ns'));
if size(d,1)>size(d,2)
    d=d';
    s=s';
end
ds = [d;s]'; %d = distance, s=correlation
ds(isnan(ds(:,2)),:)=[];
clear_all_but('ds','a2nd','a2ns');
optM = [20,20];
[temphist,bins]=CalcNormHistByDistances(ds,a2nd,a2ns,optM);
counts = repmat(sum(temphist,2),[1,optM(1)]);
counts = counts/sum(counts(:,1));
counts2 = repmat(sum(temphist,1),[optM(2),1]);
counts2 = counts2/sum(counts2(1,:));
expected = counts.*counts2;
expected(expected==0)=inf;
normhist = temphist./expected;
%------Plot--------%
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
s = surface(xo,yo,normhist);
az=0;%35
el=90;%32
view([az,el]);
set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
axis tight;
grid off;
hCbar = colorbar(gca);
set(gca,'PlotBoxAspectRatio',[1,1,1],'TickDir','out','FontSize',9);
ylabel('maximum astrocyte-neuron correlation');
xlabel('distance [um]');
ylabel(hCbar,'probability');
set(hCbar,'TickDirection','out');

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
clear_all_but('dl','a2nd','a2nL');
optM = [20,20];

[temphist,bins]=CalcNormHistByDistances(dl,a2nd,a2nL,optM);
counts = repmat(sum(temphist,2),[1,optM(1)]);
counts = counts/sum(counts(:,1));
counts2 = repmat(sum(temphist,1),[optM(2),1]);
counts2 = counts2/sum(counts2(1,:));
expected = counts.*counts2;
expected(expected==0)=inf;
normhist = temphist./expected;
normhist = normhist/trapz(bins{2},trapz(bins{1},normhist));
%---------Plot---------%
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
s = surface(xo,yo,normhist);
az=0;%35
el=90;%32
view([az,el]);
set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
axis tight;
grid off;
hCbar = colorbar(gca);
caxis([nanmin(normhist(:)),nanmax(normhist(:))]);
set(gca,'PlotBoxAspectRatio',[1,1,1],'TickDir','out','FontSize',9);
ylabel('astrocyte-neuron lag at maximum correlation [ms]');
xlabel('distance [um]');
ylabel(hCbar,'incidence [counts]');
set(hCbar,'TickDirection','out');
