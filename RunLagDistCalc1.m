mf = matfile('F:\CosSimTemp2.mat');
MaxCosSim = cell(9,1);

for i=1:9
%---Max Corr and Lags (at max corr)---%
    temp = mf.CosSim(1,i);
    temp = temp{1};
    [MaxCosSim{i},ix{i}] = nanmax(temp,[],3);
%     MaxCosSim{i}=ms;
end
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')
clear temp; 
LagAtMax = cell(9,1);
parfor i=1:9
%     LagAtMax{i} = tril(nan(size(MaxCosSim{i})))+triu(lags{i}(ix{i}),1);
%---Calc lag vs. distance---%
%     nnL = LagAtMax{i}(1:size(DataSet{i}.ic,2),1:size(DataSet{i}.ic,2));
%     aaL = LagAtMax{i}(size(DataSet{i}.ic,2)+1:end,size(DataSet{i}.ic,2)+1:end);
    anL = LagAtMax{i}(1:size(DataSet{i}.ic,2),size(DataSet{i}.ic,2)+1:end);
    [n2nd{i},n2nL{i}] = CalcValueDist(nnL,DataSet{i}.ic,[],[],'n2n');
    [a2ad{i},a2aL{i}] = CalcValueDist(aaL,[],DataSet{i}.mask,[],'a2a');
    [a2nd{i},a2nL{i}] = CalcValueDist(anL,DataSet{i}.ic,DataSet{i}.mask,DataSet{i}.channel,'a2n');
%     a2a = [a2ad{i}(~isnan(a2aL{i})),a2aL{i}(~isnan(a2aL{i}))]';
%     a2n = [a2nd{i}(~isnan(a2nL{i})),a2nL{i}(~isnan(a2nL{i}))]';
%     n2n = [n2nd{i}(~isnan(n2nL{i})),n2nL{i}(~isnan(n2nL{i}))]';
    %     a2a = a2a(:,a2a(2,:)~=0);
    %     a2aB{i} = nsOPTBINS(a2a);
    %     a2n = a2n(:,a2n(2,:)~=0);
    %     a2nB{i} = nsOPTBINS(a2n);
    %     n2n = n2n(:,n2n(2,:)~=0);
    %     n2nB{i} = nsOPTBINS(n2n);
end
save('LagDist2.mat','a2ad','a2aL','n2nd','n2nL','a2nd','a2nL','-v7.3');%,'a2aB','a2nB','n2nB');
save('CorrDistributions2.mat','LagAtMax','MaxCosSim','-v7.3');

%---Calc Corr vs. Distance---%
clear_all_but('DataSet','MaxCosSim');
parfor i=1:9
    t= DataSet{i}.t;
    ic= DataSet{i}.ic;
    mask = DataSet{i}.mask;
    a2a=MaxCosSim{i}(size(ic,2)+1:end,size(ic,2)+1:end);
    n2n=MaxCosSim{i}(1:size(ic,2),1:size(ic,2));
    a2n=MaxCosSim{i}(1:size(ic,2),size(ic,2)+1:end);
    ch = DataSet{i}.channel;
    [a2nd{i},a2ns{i}] = CalcValueDist(a2n,ic,mask,ch,'a2n');
    [a2ad{i},a2as{i}] = CalcValueDist(a2a,[],mask,[],'a2a');
    [n2nd{i},n2ns{i}] = CalcValueDist(n2n,ic,[],[],'n2n');
    % a2aB{i} = nsOPTBINS([a2ad{i}(~isnan(a2as{i})),a2as{i}(~isnan(a2as{i}))]');
    % a2nB{i} = nsOPTBINS([a2nd{i}(~isnan(a2ns{i})),a2ns{i}(~isnan(a2ns{i}))]');
    % n2nB{i} = nsOPTBINS([n2nd{i}(~isnan(n2ns{i})),n2ns{i}(~isnan(n2ns{i}))]');
end
save('CorrDist2.mat','a2nd','a2ns','a2ad','a2as','n2nd','n2ns');%,'a2aB','a2nB','n2nB');