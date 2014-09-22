function factor = Calc2DhistNormalization(ic,binsd,MeaMap,mask,which)
%Function which takes the map of electrode positions(MeaMap), the list of
%channels in the recording (ic) and the bins with which to calculate the
%histogram. Outputs the fraction of distances measurable for each bin and
%for each electrode.
nnB = size(binsd,2);
switch which
    case 'n2n'
        DistMat = CreateElectrodeDistanceTable();
        factor=zeros(numel(ic),nnB(1));
        for j=1:numel(ic)
            for k=1:nnB(1)-1
                factor(j,k) = NormalizeDistanceCounts(DistMat,MeaMap,ic(j),binsd(k),binsd(k+1));
            end
            factor(j,:) = factor(j,:)./trapz(binsd,factor(j,:));    
        end

        combs = VChooseK(1:numel(ic),2);
        factor = factor(combs(:,1),:);
        %--Normalization--%
%         factor = factor./repmat(nansum(factor,1),size(factor,1),1);
        % This means a 1 is when every possible pair has a tick in a given
        % bin. 
        factor(isnan(factor))=0;
    case 'a2n'
        DistMat = CreateElectrodeDistanceTable();
        factor=zeros(1,nnB(1));
        for k=1:nnB(1)-1
            factor(k) = NormalizeDistanceCounts(DistMat,MeaMap,ic,binsd(k),binsd(k+1));
        end
        factor = repmat(factor./trapz(binsd,factor),mask,1);
%         factor = factor./(255*size(factor,1));
%         factor = factor./repmat(nansum(factor,1),size(factor,1),1);
        factor(isnan(factor))=0;
    case 'a2a'
        distmat = triu(regiondist(mask),1);
        factor=zeros(size(distmat,1),nnB(1));
        for i=1:size(distmat,1)
            for k=1:nnB-1
                factor(i,k) =  sum((distmat(i,:)>=binsd(k))&(distmat(i,:)<binsd(k+1)));
            end
        end
        combs = VChooseK(1:size(distmat,1),2);
        factor = factor(combs(:,1),:);
        factor = factor./repmat(nansum(factor,1),size(factor,1),1);
        factor(isnan(factor))=0;
    otherwise
        error('Not an appropriate selection. The options are: (1) A2N, (2) N2N or (3) A2A');
end

end