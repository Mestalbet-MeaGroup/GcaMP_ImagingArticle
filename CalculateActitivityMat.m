function [ActivityMat,orderedAmat]=CalculateActitivityMat(MaxCosSim,traces,fr,ic,nx,ax)
%Function which calculates a matrix arranged like MaxCosSim except instead
%of correlation values, shows the geometric average activity (number of
%spikes or astro peaks). Also returns the matrix dendrogram ordered by DendrogramOrderMaxCosSim

ActivityMat = zeros(size(MaxCosSim));
[~,pp,~]=CalcPeakStartEnd(traces); 
for i=1:size(ic,2)
    for j=1:size(traces,2)
        ActivityMat(i,j+size(ic,2)) = sqrt(nanmax(fr(:,i))*numel(pp{j}));
    end
end
ActivityMat=ActivityMat+ActivityMat';

for i=1:size(ic,2)
    for j=1:size(ic,2)
        ActivityMat(i,j) = sqrt(nanmax(fr(:,i))*nanmax(fr(:,j)));
    end
end

for i=1:size(traces,2)
    for j=1:size(traces,2)
        ActivityMat(i+size(ic,2),j+size(ic,2)) = sqrt(numel(pp{i})*numel(pp{j}));
    end
end


a2a=ActivityMat(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=ActivityMat(1:size(ic,2),1:size(ic,2));
a2n=ActivityMat(1:size(ic,2),size(ic,2)+1:end);

a2a = a2a+a2a'+eye(size(a2a));
n2n = n2n+n2n'+eye(size(n2n));

orderedAmat = zeros(size(MaxCosSim));
orderedAmat(1:numel(nx),1:numel(nx))=n2n(nx,nx);
orderedAmat(numel(nx)+1:numel(nx)+numel(ax),numel(nx)+1:numel(nx)+numel(ax))=a2a(ax,ax);
orderedAmat(1:numel(nx),numel(nx)+1:numel(nx)+numel(ax))=a2n(nx,ax);
orderedAmat(numel(nx)+1:numel(nx)+numel(ax),1:numel(nx))=a2n(nx,ax)';
end