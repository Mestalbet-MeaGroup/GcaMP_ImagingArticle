load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
parfor i=1:9
    t= DataSet{i}.t;
    ic= DataSet{i}.ic;
    traces = DataSet{i}.dfTraces;
    time=DataSet{i}.dfTime;
    bs = DataSet{i}.bs;
    be = DataSet{i}.be;
%     [atrigpsth{i},atriglag{i}]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);
    [ntrigpsth{i},ntriglag{i}]=CalcPSTHneurotriggers(t,ic,traces,time,bs,be);
    %     clear_all_but('atrigpsth','atriglag','ntrigpsth','ntriglag','i');
end
% clear_all_but('atrigpsth','atriglag');
% save('PSTHdata_byStarts.mat','atrigpsth','atriglag');
% save('PSTHdata_byPeaks.mat','atrigpsth','atriglag');
save('PSTHdata_byPeaks.mat','ntrigpsth','ntriglag','-append');








%
% load('PSTHdata_byPeaks.mat')
subplot(2,1,1)
hold all;
atrigmean =[];
for i=1:9
   
    plot(atriglag{i},nansum(zscore(atrigpsth{i},0,2),1));
    atrigmean = [atrigmean;nansum(zscore(atrigpsth{i},0,2))];
end
plot(atriglag{i},nanmean(atrigmean,1),'-k','LineWidth',5);
title('Astrocyte Triggered PSTH');
subplot(2,1,2)
hold all;
for i=1:9
   plot(ntriglag{i},nansum(zscore(ntrigpsth{i},0,2),1));
end
title('Burst Triggered Astrocyte Peak PSTH');