function [n2nL,a2aL,a2nL,ic]=CalcLagsAtMaxCorr(i,ic)
% i=1;
% CS  = matfile('CosSimTemp2.mat');
% db = matfile('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
load('F:\CosSimTemp2.mat','lags');
load('CorrDistributions2.mat', 'LagAtMax')
%---fs---%
% fs = DataSet{i,1}.fs;
% fs = fs{1}.fs;
%--------%

%---cs---%
% csi = CS.CosSim(1,i);
% csi = csi{1};
%--------%

%---ic---%
% ic = DataSet{i,1}.ic;
% ic = ic{1}.ic;
%--------%
% clear CS; clear DataSet;
%---Lag Mat---%
% lags = (-size(csi,3)/2:size(csi,3)/2-1).*(1/fs);
% [csmax,idx] = nanmax(csi,[],3);
%-------------%
lags = tril(nan(size(LagAtMax{i})))+triu(lags{i}(LagAtMax{i}),1);

n2nL = lags(1:size(ic,2),1:size(ic,2));
a2aL = lags(size(ic,2)+1:end,size(ic,2)+1:end);
a2nL = lags(1:size(ic,2),size(ic,2)+1:end);
end