% Calculate Astrocyte to Neural correlation by using burst start times.
fclose('all');clear all; close all;clc;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat')
%% Astrocyte transients time versus binary burst starts
% maxlag=5000;
% for i=1:size(DataSet,1);
%     if isfield(DataSet{i},'sb_bs')
%         bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs])./12000;
%         be = sort([DataSet{i}.be,DataSet{i}.sb_be])./12000;
%     else
%         bs = DataSet{i}.bs./12000;
%         be = DataSet{i}.be./12000;
%     end
%     bs(bs<=DataSet{i}.dfTime(1))=[];
%     bs(bs>=DataSet{i}.dfTime(end))=[];
%     [time,ix]=sort([DataSet{i}.dfTime,bs]);
%     for j=1:size(DataSet{i}.dfTraces,2)
%         zval = interp1(DataSet{i}.dfTime,zscore(DataSet{i}.dfTraces(:,j)),bs,'spline');
%         trace = [zscore(DataSet{i}.dfTraces(:,j))',zval];
%         trace = trace(ix);
%         bursts = zeros(size(trace));
%         bursts(find(ismember(time,bs)))=1;
%         Corr{i}(j,:) = xcorr(trace,bursts,maxlag,'coeff');
%     end
% end
% lags = -maxlag:maxlag;
% 
%% Astrocyte transients versus burst start with height determined by FR of burst
% maxlag=5000;
% for i=1:size(DataSet,1);
%     numc=[];
%     if isfield(DataSet{i},'sb_bs')
%         bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs])./12000;
%         be = sort([DataSet{i}.be,DataSet{i}.sb_be])./12000;
%     else
%         bs = DataSet{i}.bs./12000;
%         be = DataSet{i}.be./12000;
%     end
%     be(bs<=DataSet{i}.dfTime(1) | bs>=DataSet{i}.dfTime(end))=[];
%     bs(bs<=DataSet{i}.dfTime(1) | bs>=DataSet{i}.dfTime(end))=[];
%     [time,ix]=sort([DataSet{i}.dfTime,bs]);
%     [t,it]=sort(DataSet{i}.t./12000);
%     vec=ConvertIC2Samora(DataSet{i}.ic);
%     vec=sparse(vec(it));
%     mat=bsxfun(@gt,sparse(t),sparse(bs'))&bsxfun(@lt,sparse(t),sparse(be')); %number of spikes in a burst
%     vec = repmat(vec,[size(mat,1),1]);
%     channels = vec.*mat;
%     for k=1:size(channels,1), numc(k) = numel(unique(channels(k,:))); end;
%     fr = (sum(mat,2)./(be-bs)')./numc';
%     bursts = zeros(size(time));
%     bursts(ismember(time,bs))=fr;
%     for j=1:size(DataSet{i}.dfTraces,2)
%         zval = interp1(DataSet{i}.dfTime,zscore(DataSet{i}.dfTraces(:,j)),bs,'spline');
%         trace = [zscore(DataSet{i}.dfTraces(:,j))',zval];
%         trace = trace(ix);
%         CorrFR{i}(j,:) = xcorr(trace,zscore(bursts),maxlag,'coeff');
%     end
% end
% lags = -maxlag:maxlag;

%---Load instead of calculate---%
load('CorrFR.mat');

%% Plot Differences
% for i=1:size(DataSet,1)
%     RemainingCorr{i} = CorrFR{i}-Corr{i};
% end
% i=8;
% subplot(3,1,1)
% imagesc(lags,1:size(CorrFR{i},2),CorrFR{i}); colorbar;
% title('With FR');
% 
% subplot(3,1,2)
% imagesc(lags,1:size(Corr{i},2),Corr{i}); colorbar;
% title('Binary');
% 
% subplot(3,1,3)
% imagesc(lags,1:size(RemainingCorr{i},2),RemainingCorr{i}); colorbar;
% title('Difference');
% 
%% Pairwise correlation between astrocyte transients vs burst start with height determined by FR of each electrode
maxlag=5000;
CorrFR=cell(9,1);
for i=1:size(DataSet,1);
    if isfield(DataSet{i},'sb_bs')
        bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs])./12000;
        be = sort([DataSet{i}.be,DataSet{i}.sb_be])./12000;
    else
        bs = DataSet{i}.bs./12000;
        be = DataSet{i}.be./12000;
    end
    be(bs<=DataSet{i}.dfTime(1) | bs>=DataSet{i}.dfTime(end))=[];
    bs(bs<=DataSet{i}.dfTime(1) | bs>=DataSet{i}.dfTime(end))=[];
    [time,ix]=sort([DataSet{i}.dfTime,bs]);
    [t,it]=sort(DataSet{i}.t./12000);
    ic = DataSet{i}.ic;
    vec=ConvertIC2Samora(ic);
    vec=sparse(vec(it));
    mat=bsxfun(@gt,sparse(t),sparse(bs'))&bsxfun(@lt,sparse(t),sparse(be')); %number of spikes in a burst
    vec = repmat(vec,[size(mat,1),1]);
    channels = vec.*mat;
    numc=numel(unique(vec));
    bursts = zeros(numc,length(time));
    for k=1:numc
            fr = sum(channels==ic(1,k),2)./(be-bs)';
            bursts(k,ismember(time,bs))=fr;          
    end
    combs = allcomb(1:numc,1:size(DataSet{i}.dfTraces,2));
    e1 = combs(:,1);
    e2 = combs(:,2);
    temp=zeros(size(combs,1),maxlag*2+1);
    temp1=temp;
    linidx= sub2ind([numc,size(DataSet{i}.dfTraces,2)],e1,e2);
    bursts = bursts(e1,:);
    tr = DataSet{i}.dfTraces(:,e2);
    timetemp=DataSet{i}.dfTime;
    parfor j=1:size(combs,1)
        zval = interp1(timetemp,zscore(tr(:,j)),bs,'spline');
        trace = [zscore(tr(:,j))',zval];
        trace = trace(ix);
        temp(j,:) = xcorr(trace,zscore(bursts(j,:)),maxlag,'coeff');
    end
   temp1(linidx,:)=temp;
   PairWiseCorrFR{i}=reshape(temp1,[numc,size(DataSet{i}.dfTraces,2),maxlag*2+1]);
end
lags = -maxlag:maxlag;
