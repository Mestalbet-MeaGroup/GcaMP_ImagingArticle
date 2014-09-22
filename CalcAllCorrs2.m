load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
maxlag=5000;
for which=6
    r = [DataSet{which}.FR,DataSet{which}.dfTraces]';
    CosSim = CalcCorrFast(r,maxlag);
    lags  = (-maxlag:maxlag)./DataSet{which}.fs;
end
clear_all_but('CosSim','lags');
save('CosSimTemp2.mat','CosSim','lags');

[MaxCosSim,ix] = nanmax(CosSim,[],3);
LagAtMax = tril(nan(size(MaxCosSim)))+triu(lags(ix),1);