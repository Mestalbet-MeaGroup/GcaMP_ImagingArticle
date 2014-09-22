function p = PeakFinder_SW(tr)
%A wrapper for peakfinder which splits a trace into windows.

% trace = tr;
% tr = smooth(zscore(tr(:,1)),50,'sgolay',5);
xvals = 1:length(tr);
yvals = tr;

numwins=4;
nwind = round(length(tr)/numwins);
noverlap=100;

[xvals,index]=sort(xvals);
yvals=yvals(index);

nx = length(xvals);                            %# length of sequence
ncol = fix((nx-noverlap)/(nwind-noverlap));    %# number of sliding windows
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);  %# starting index of each
idx = bsxfun(@plus, (1:nwind)', colindex)-1;   %# each column is the indexes of one window. Each row is a sliding step forward.
newy=yvals(idx);

p=[];
for i=1:numwins
    temp = newy(:,i);
    [~,s1,~] = deleteoutliers(temp,0.01);
    s2 = peakfinder(temp,(max(temp)-min(temp))/10);
    idx = intersect(s1,s2);
    p = [idx+colindex(i);p];
end
p=unique(p);

around = 50;
trace = padarray(tr,[around,0],nan,'both');
p = arrayfun(@(x) find(trace==nanmax(trace(x-around:x+around)),1,'First'),p+around)-around;
p=unique(p);

end