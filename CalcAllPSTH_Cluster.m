function [atrigpsth, atriglag] = CalcAllPSTH_Cluster()
load('/home/nl1001/MdcsDataLocation/freiburg/R2013a/remote/DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');

parfor i=1:9
    t= DataSet{i}.t;
    ic= DataSet{i}.ic;
    traces = DataSet{i}.dfTraces;
    time=DataSet{i}.dfTime;
    bs = DataSet{i}.bs;
    be = DataSet{i}.be;
    [atrigpsth{i},atriglag{i}]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);
    %     [ntrigpsth{i},ntriglag{i}]=CalcPSTHneurotriggers(t,ic,traces,time,bs,be);
    %     clear_all_but('atrigpsth','atriglag','ntrigpsth','ntriglag','i');
end
end
% clear_all_but('atrigpsth','atriglag');
% save('PSTHdata_byStarts.mat','atrigpsth','atriglag');

%
% load('PSTHdata_byPeaks.mat')
% subplot(2,1,1)
% hold all;
% for i=1:9
%    plot(atriglag{i},nansum(zscore(atrigpsth{i},0,2),1));
% end
% title('Astrocyte Triggered PSTH');
% subplot(2,1,2)
% hold all;
% for i=1:9
%    plot(ntriglag{i},nansum(zscore(ntrigpsth{i},0,2),1));
% end
% title('Burst Triggered Astrocyte Peak PSTH');