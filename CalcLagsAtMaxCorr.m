function [n2nL,a2aL,a2nL,lagmat,ic]=CalcLagsAtMaxCorr(i)
% i=1;
CS  = matfile('F:\CosSimTemp.mat');
db = matfile('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');

%---fs---%
fs = db.DataSet(1,i);
fs = fs{i}.fs;
%--------%

%---cs---%
csi = CS.CosSim(i,i);
csi = csi{1};
%--------%

%---ic---%
ic = db.DataSet(i,i);
ic = ic{i}.ic;
%--------%
clear CS; clear db;
%---Lag Mat---%
lags = (-size(csi,3)/2:size(csi,3)/2-1).*(1/fs);
[csmax,idx] = nanmax(csi,[],3);
%-------------%

lagmat = lags(idx);
n2nL = lagmat(1:size(ic,2),1:size(ic,2));
a2aL = lagmat(size(ic,2)+1:end,size(ic,2)+1:end);
a2nL = lagmat(1:size(ic,2),size(ic,2)+1:end);
end