function DistMat = CreateElectrodeDistanceTable();
MeaImage = [];
MeaMap=[];
load('MeaMapPlot.mat','MeaImage','MeaMap');
MeaMap = rot90(MeaMap,-1);
overlay = otsu(imresize(MeaImage,1.0441),45);
overlay(overlay>1)=0;
overlay=bwareaopen(overlay, 500);
overlay = (overlay>=1);
%--Erase manually two spots---%
overlay(1750:1830,2705:2770)=0;
overlay(1600:1720,2745:2830)=0;
%-----------------------------%
centers = imfindcircles(overlay,[12,35],'ObjectPolarity','bright','method','TwoStage','Sensitivity',0.70);
ul = [132.984326492796,3104.01444974662];
ur = [3132.77588690862,3088.04457677175];
ll = [117.990815075410,105.047826670482];
lr = [3119.92295829176,89.3609937082245];
centers = cat(1,centers,ul,ur,ll,lr);
centers=sortrows(centers,2);

% [xpos,a]=sort(centers(:,2));
% ypos=centers(a,1);
xpos = centers(:,2);
ypos = centers(:,1);
xpos=reshape(xpos,16,16);
ypos=reshape(ypos,16,16);

DistMat = squareform(pdist(centers));

%----How to look up a distance between two electrodes in the lookup table, DistMat---%
% [x1,y1] = find(MeaMap==ic(source));
% [x2,y2] = find(MeaMap==ic(target));
% ind1 = sub2ind(size(MeaMap),x1,y1);
% ind2 = sub2ind(size(MeaMap),x2,y2);
% distance1_2 = DistMat(ind1,ind2);


end