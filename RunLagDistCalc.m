load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat')
parfor i=1:9
    [nnL,aaL,anL,~]=CalcLagsAtMaxCorr(i,DataSet{i}.ic);
    [n2nd{i},n2nL{i}] = CalcValueDist(nnL,DataSet{i}.ic,[],[],'n2n');
    [a2ad{i},a2aL{i}] = CalcValueDist(aaL,[],DataSet{i}.mask,[],'a2a');
    [a2nd{i},a2nL{i}] = CalcValueDist(anL,DataSet{i}.ic,DataSet{i}.mask,DataSet{i}.channel,'a2n');
    a2a = [a2ad{i}(~isnan(a2aL{i})),a2aL{i}(~isnan(a2aL{i}))]';
    a2aB{i} = nsOPTBINS(a2a);
    a2n = [a2nd{i}(~isnan(a2nL{i})),a2nL{i}(~isnan(a2nL{i}))]'
    a2nB{i} = nsOPTBINS(a2n);
    n2n = [n2nd{i}(~isnan(n2nL{i})),n2nL{i}(~isnan(n2nL{i}))]'
    n2nB{i} = nsOPTBINS(n2n);
end
save('LagDist.mat','a2ad','a2aL','n2nd','n2nL','a2nd','a2nL','a2aB','a2nB','n2nB');

% clear all;
% load('LagDist.mat','n2nL','a2aL','a2nL');
% parfor i=1:9
% [aaHisL{i},aaBinL{i}]=CalcHist(a2aL{i});
% [anHisL{i},anBinL{i}]=CalcHist(a2nL{i});
% [nnHisL{i},nnBinL{i}]=CalcHist(n2nL{i});
% end
% save('LagHistos.mat');

% clear all; 
% for i=1:9
% loadcult(i);
% [a2nd{i},a2ns{i}] = CalcValueDist(a2n,ic,mask,ch,'a2n');
% [a2ad{i},a2as{i}] = CalcValueDist(a2a,[],mask,[],'a2a');
% [n2nd{i},n2ns{i}] = CalcValueDist(n2n,ic,[],[],'n2n');
% a2aB{i} = nsOPTBINS([a2ad{i}(~isnan(a2as{i})),a2as{i}(~isnan(a2as{i}))]');
% a2nB{i} = nsOPTBINS([a2nd{i}(~isnan(a2ns{i})),a2ns{i}(~isnan(a2ns{i}))]');
% n2nB{i} = nsOPTBINS([n2nd{i}(~isnan(n2ns{i})),n2ns{i}(~isnan(n2ns{i}))]');
% clear_all_but('a2nd','a2ns','a2ad','a2as','n2nd','n2ns','a2aB','a2nB','n2nB','i');
% end
% save('CorrDist.mat','a2ad','a2as','n2nd','n2ns','a2nd','a2ns','a2aB','a2nB','n2nB');