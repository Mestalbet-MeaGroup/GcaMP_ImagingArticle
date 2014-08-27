function d = FindMinDifferences(ts1,ts2);
% Finds the minimum distance between every element of ts1 and every element of ts2
 m = size(ts1,2); n = size(ts2,2);
 [~,p] = sort([ts1,ts2]);
 q = 1:m+n; q(p) = q;
 t = cumsum(p>m);
 r = 1:n; r(t(q(m+1:m+n))) = r;
 s = t(q(1:m));
 id = r(max(s,1));
 iu = r(min(s+1,n));
 [d,~] = min([abs(ts1-ts2(id));abs(ts2(iu)-ts1)]);
%  ib = id+(it-1).*(iu-id);
end