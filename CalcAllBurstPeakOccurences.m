load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')
parfor k=1:size(DataSet,1)
    [numBI(k),numNoB(k),numSB(k),PeaksWithBursts{k},BurstsWithPeaks{k},PeaksNoBursts{k},BurstsNoPeaks{k},PeaksWithSBs{k},SBsWithPeaks{k},peakIndex{k}]=CalcBurstPeakOccurenceRate(DataSet,k);
end
save('BurstPeakOccurance.mat');