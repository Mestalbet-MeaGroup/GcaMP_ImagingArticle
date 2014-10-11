function [peakIndex,OnIndex,OffIndex,ClosestValues,BurstIndex]=FindPeaksNearestBursts(traces,time,BurstOnsets,t);
% Function which take traces, times and burst starts and (1) Finds the
% peaks within each trace and (2) Calculates the minimum distance between
% each peak and its nearest burst start.  

[OnIndex,peakIndex,OffIndex]=CalcPeakStartEnd(traces); 
 b = arrayfun(@(x) max(t(find(t>=x,10,'First'))),BurstOnsets)./12000; %convert to seconds
for j=1:size(traces,2)
    
    a = time(peakIndex{j}); %time of peak in seconds
    m = numel(a);
    n = numel(b);
    
    [~,p] = sort([a,b]);
    q = 1:m+n;
    q(p) = q;
    t = cumsum(p>m);
    r = 1:n; r(t(q(m+1:m+n))) = r;
    s = t(q(1:m));
    id = r(max(s,1));
    iu = r(min(s+1,n));
    ClosestValues{j} = (a-b(id));  
    BurstIndex{j} = id;
end


end