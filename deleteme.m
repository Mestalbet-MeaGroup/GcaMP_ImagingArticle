parfor i=1:9
    if ~isempty(DataSet{i}.sbs)
        DataSet{i}.sb_bs=[];
        DataSet{i}.sb_be=[];
        for j=1:numel(DataSet{i}.sbs)
            [bs,be] = FindBurstsWithinSBs_special(DataSet{i}.t,DataSet{i}.ic,DataSet{i}.sbs(j),DataSet{i}.sbe(j),[],[],[],0);
            DataSet{i}.sb_bs=[DataSet{i}.sb_bs,bs-DataSet{i}.sbs(j)];
            DataSet{i}.sb_be=[DataSet{i}.sb_be,be-DataSet{i}.sbs(j)];
        end
    end
end

parfor k=1:9
[numBI{k},numNoB{k},numSB{k},PeaksWithBursts{k},BurstsWithPeaks{k},PeaksNoBursts{k},BurstsNoPeaks{k},PeaksWithSBs{k},SBsWithPeaks{k},peakIndex{k}]=CalcBurstPeakOccurenceRate(DataSet,k);
end
clear DataSet;
save('BurstPeakOccurance_withSBbursts_onlyforward.mat');

plotyy(bins{1},mean(test(1:end/2,:),1),bins{1},mean(test(end/2+1:end,:),1));
plot(bins{1},mean(test(end-end/3:end,:),1)./mean(test(1:end/3,:),1),'.-');