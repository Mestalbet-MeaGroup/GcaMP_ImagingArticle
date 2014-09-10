function factor = Calc2DhistNormalization(ic,bins,MeaMap,mask,which)
%Function which takes the map of electrode positions(MeaMap), the list of
%channels in the recording (ic) and the bins with which to calculate the
%histogram. Outputs the fraction of distances measurable for each bin and
%for each electrode.
nnB = size(bins,2);
switch which
    case 'n2n'
        DistMat = CreateElectrodeDistanceTable();
        factor=zeros(numel(ic),nnB(1));
        for j=1:numel(ic)
            for k=1:nnB(1)-1
                factor(j,k) = NormalizeDistanceCounts(DistMat,MeaMap,ic(j),bins(k),bins(k+1));
            end
        end
        factor(isinf(factor))=0;
%         factor = factor./sum(factor,2);
        combs = VChooseK(1:numel(ic),2);
        factor = factor(combs(:,1),:);
    case 'a2n'
        DistMat = CreateElectrodeDistanceTable();
        factor=zeros(1,nnB(1));
        for k=1:nnB(1)-1
            factor(k) = NormalizeDistanceCounts(DistMat,MeaMap,ic,bins(k),bins(k+1));
        end
        factor = repmat(factor,mask,1);
    case 'a2a'
        distmat = triu(regiondist(mask),1);
        factor=zeros(size(distmat,1),nnB(1));
        for i=1:size(distmat,1)
            for k=i:nnB-1
                factor(i,k) =  1/sum((distmat(i,:)>=bins(k))&(distmat(i,:)<bins(k+1)));
            end
        end
        factor(isinf(factor))=0;
        combs = VChooseK(1:size(distmat,1),2);
        factor = factor(combs(:,1),:);
    otherwise
        error('Not an appropriate selection. The options are: (1) A2N, (2) N2N or (3) A2A');
end

end