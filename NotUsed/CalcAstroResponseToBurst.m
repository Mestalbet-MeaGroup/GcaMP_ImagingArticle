function percent=CalcAstroResponseToBurst(DataSet,k)
% Calculates the percentage of peaks that are within 1 second of a burst
% start out of the sum of peaks within 1 second and greater than 3 seconds

[~,peakIndex,~]=CalcPeakStartEnd(DataSet{k}.dfTraces);

%%
counter=1;cf=1;
for i=1:size(peakIndex,2)
    peaks = DataSet{k}.dfTime(peakIndex{i});
    distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
    closeby = zeros(size(distances));
    closeby = abs(distances)<1;
    far = zeros(size(distances));
    far = abs(distances)<=3;
    [pc{i},~]=find(closeby==1);
    pf{i}=find(sum(far==0,2)==size(distances,2));
    
    percent(i) = numel(pc{i}) / (numel(pc{i})+numel(pf{i}));
end
end