load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
t=DataSet{6}.ic;
ic=DataSet{6}.ic;
t=DataSet{6}.t;
bs = DataSet{6}.bs;
be = DataSet{6}.be;
eventS = bs(find((be-bs)==max(be-bs)));
eventE = be(find((be-bs)==max(be-bs)));
[t,ix]=sort(t,'ascend');
vec=ConvertIC2Samora(ic);
vec=vec(ix);
vec = vec((t>=eventS)&(t<=eventE));
e1 = unique(vec,'stable');
eventFR = arrayfun(@(x) sum(vec==x), 1:max(e1));
PlotBurstPropogationAnimation(e1,eventFR);


%-------------%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
t = DataSet{6}.t;
ic = DataSet{6}.ic;
sbs = DataSet{6}.sbs;
sbe = DataSet{6}.sbe;
eventS= sbs(2);
eventE= sbe(2);
[t,ix]=sort(t,'ascend');
vec=ConvertIC2Samora(ic);
vec=vec(ix);
vec = vec((t>=eventS)&(t<=eventE));
e1 = unique(vec,'stable');
eventFR = arrayfun(@(x) sum(vec==x), 1:max(e1));
PlotBurstPropogationAnimation(e1,eventFR);


%------------%
fclose('all');clear all; close all;clc;
load('DataSetOpto_trim4HA.mat');
k=1; kk=2; 

t =  DataSetStims{k}.t;
ic = DataSetStims{k}.ic;
bs = DataSetStims{k}.sbs;
be = DataSetStims{k}.sbe;
eventS = bs(find((be-bs)==max(be-bs)));
eventE = be(find((be-bs)==max(be-bs)));

% 
% [bs1,be1] = FindBurstsWithinSBs(DataSetStims{k}.Trim.t,DataSetStims{k}.Trim.ic,DataSetStims{k}.sbs(kk),DataSetStims{k}.sbe(kk),k,kk,0,0);
% eventS = bs1(4);
% eventE = be1(4);

[t,ix]=sort(t,'ascend');
vec=ConvertIC2Samora(ic);
vec=vec(ix);
vec = vec((t>=eventS)&(t<=eventE));
e1 = unique(vec,'stable');
eventFR = ones(max(e1),1);
PlotBurstPropogationAnimation(e1,eventFR,10);


% [X,Y] = meshgrid(1:1:16);
% val=linspace(1,255,4);
% Z =[ones(8,8)*val(1),ones(8,8)*val(2);ones(8,8)*val(3),ones(8,8)*val(4)];
% Z1 =  filtfilt(MakeGaussian(0,2,3),1,Z);
% Z1 =  filtfilt(MakeGaussian(0,2,3),1,Z1');
% numVals = numel(unique(Z1(:)));
% 
% cmap = colormap(cbrewer('div', 'Spectral', numVals));