k=6;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat');
load('PSTHdata_byPeaks.mat');
load('BurstPeakOccurance_withSBbursts_onlyforward.mat');
pi  = peakIndex{k};
pb  = PeaksWithBursts{k};
c=1;
for i=1:size(pi,2)
    tr  = DataSet{k}.dfTraces(:,i)-min(DataSet{k}.dfTraces(:,i));
    if ~isempty(pb{i})
        pb{i}=unique(pb{i});
        for j=1:numel(pb{i})
           CutoutNear(c,:) = CreateTraceCutOuts(tr,pi{i}(pb{i}(j)),50,50);
           c=c+1;
        end
    end
end
numlevs=4;
ampsC  = nanmax(CutoutNear,[],2);
ranksC = otsu(ampsC,numlevs);
CreatePSTHFigure(DataSet,atriglag,atrigpsth,6,CutoutNear,ranksC)