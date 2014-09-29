function order = SortPSTHpairs(psth);

% [PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);

% clusters = kmeans(PairWisePSTH',5);
% [~,ix]=max([sum(clusters==1),sum(clusters==2),sum(clusters==3),sum(clusters==4),sum(clusters==5)]);
% clusters(clusters~=ix)=0;
% clusters(clusters==ix)=1;
% [sIDX,eIDX]=initfin(~clusters');
% for i=1:numel(sIDX)
%     sc(i) = eIDX(i)-sIDX(i);
% end
% [~,ix]=nanmax(sc);
% 
% psth = PairWisePSTH(:,sIDX(ix):eIDX(ix));

for j=1:size(psth,1)
    [corr,lag]=xcorr(zscore(psth(1,:)),zscore(psth(j,:)),100,'unbiased');
    offset(j) = lag(find(corr==max(corr),1,'First'));
end
[~,order]=sort(offset);

end