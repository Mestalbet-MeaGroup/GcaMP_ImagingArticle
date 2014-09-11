function factor = NormalizeDistanceCounts(DistMat,MeaMap,elec,bins,bine);
%Function which accepts an electrode id and looks up the number of
%distances between it and all other electrodes which are at distance, d. 
% eidx = the index number of this electrode 
% DistMat - produced from CreateElectrodeDistanceTable using the scatter
% plot of MEA locations
% bins  = previous bin
% bine = the next bin
[x1,y1] = find(MeaMap==elec);
eidx = sub2ind(size(MeaMap),x1,y1);
dists = DistMat(eidx,:);
factor =sum((dists>=bins)&(dists<=bine)); %or the inverse, not sure.
% factor =1/sum((dists>=bins)&(dists<bine)); %or the inverse, not sure.
end