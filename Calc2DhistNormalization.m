function factor = Calc2DhistNormalization(ic,bins,MeaMap)

DistMat = CreateElectrodeDistanceTable();
nnB = size(bins,2);
factor=zeros(numel(ic),nnB(1));
for j=1:numel(ic)
    for k=1:nnB(1)-1
        factor(j,k) = NormalizeDistanceCounts(DistMat,MeaMap,ic(j),bins(k),bins(k+1));
%           factor(j,k) = 1;
    end
end
combs = VChooseK(1:numel(ic),2);
factor = factor(combs(:,1),:);
% f1 = factor(combs(:,1),:);
% f2 = factor(combs(:,2),:);
% factor =f1+f2;

end