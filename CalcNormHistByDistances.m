function [normhist,bins]=CalcNormHistByDistances(pooled,bycultdist,bycultscore,optM)
load('MeaMapPlot.mat','MeaMap');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')

bin = {linspace(min(pooled(:,1)),max(pooled(:,1)),optM(1)),linspace(min(pooled(:,2)),max(pooled(:,2)),optM(2))};
normhist=zeros(numel(bin{1}),numel(bin{2}));
DistMat = CreateElectrodeDistanceTable();
for i=1:size(DataSet,1)
    [hist,bins]  = hist3([bycultdist{i}(~isnan(bycultdist{i})),bycultscore{i}(~isnan(bycultdist{i}))],'Edges',bin);
    nnB = size(bins{1},2);
    factor=zeros(1,nnB(1));
    for k=1:nnB(1)-1
        factor(k) = NormalizeDistanceCounts(DistMat,MeaMap,DataSet{i}.channel,bins{1}(k),bins{1}(k+1));
    end
    factor(factor==0)=inf;
    temphist = hist./repmat(factor,size(normhist,2),1)';
    normhist=normhist+temphist;%/trapz(bins{2},trapz(bins{1},temphist));
end
normhist = normhist./size(DataSet,1);
end