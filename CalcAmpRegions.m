function [NearFar,OnlyNear,OnlyFar]=CalcAmpRegions(NearROIs,FarROIs,r)

a=unique(NearROIs(NearROIs(:,1)==r,2));
b=unique(FarROIs(FarROIs(:,1)==r,2));
[NearFar, ia,ib]=intersect(a,b);
OnlyNear = a(setdiff(1:numel(a),ia));
OnlyFar = b(setdiff(1:numel(b),ib));
end