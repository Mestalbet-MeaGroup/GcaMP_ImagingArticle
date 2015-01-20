%Plots Pairwise astrocyte to neuron correlation versus astrocyte
%correlation to the meanfield
%% Load Data
fclose('all');clear all; close all;clc;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat')
load('CorrDistributions2.mat', 'MaxCosSim');
%---Remove Culture 4 which has virtually no transients---%
DataSet(4)=[];
MaxCosSim(4)=[];
%%
%------Create Variable Names------%
nnscores = [];
nascores = [];
astro    = [];
%---------------------------------%

for i=1:size(DataSet,1)
    corrs = MaxCosSim{i};
    corrs = triu(corrs,1)+tril(nan(size(corrs)),1);
    numc = size(DataSet{i}.ic,2);
    n2n = corrs(1:numc,1:numc);
    a2a = corrs(numc+1:end,numc+1:end);
    a2n = corrs(1:numc,numc+1:end);
    
%     nnscores = [nnscores;n2n(:)];
%     nascores = [nascores;a2n(:)];
    nascores = [nascores;nanmax(a2n,[],1)'];
    mtr= nanmean(DataSet{i}.dfTraces,2);
    temp=zeros(1,size(a2n,2));
    for j=1:size(a2n,2)
        temp(j) = CalcMaxCorr(DataSet{i}.dfTraces(:,j),mtr);
    end
    
%     temp = repmat(temp,[size(a2n,1),1]);

    astro = [astro;temp(:)];
end
remove = find(isnan(nascores));
astro(remove)=[];
nascores(remove)=[];
PlotColWiseNormHist(astro,nascores,[50,50],'both');
