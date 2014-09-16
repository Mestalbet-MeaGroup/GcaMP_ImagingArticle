function CalcAllCorrsCluster();
load('/home/nl1001/MdcsDataLocation/freiburg/R2013a/remote/DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
maxlag=8000;
for which=1
    r = [DataSet{which}.FR,DataSet{which}.dfTraces]';
    CosSim{which} = CalcCorrFast(r,maxlag);
%     lags{which} = (-maxlag:maxlag)./DataSet{which}.fs;
end
% clear_all_but('CosSim','lags');
save('/home/nl1001/MdcsDataLocation/freiburg/R2013a/remote/CosSimTemp2.mat','CosSim','-v7.3');
end