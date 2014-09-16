function [peakIndex,OnIndex,OffIndex,ClosestValues]=FindPeaksNearestBursts(traces,time,BurstOnsets);

% f = savgol(35,5,0);
[ps,peakIndex,pe]=CalcPeakStartEnd(traces);
for j=1:size(traces,2)

    bs = ps{j};
    be = pe{j};
    OnIndex{j}=bs;
    OffIndex{j}=be;
    clear bs; clear be;
    
    a = time(peakIndex{j});
    b = BurstOnsets./12000;
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
    [ClosestValues{j},~] = min([abs(a-b(id));abs(b(iu)-a)]);
end


end