function [cv,s1i,s2i]=FindNearestElements(a,b)
% Function which calculates the distance between every element of a and b
% and returns the distance between them (cv), and the index of each element
%  with the corresponding distance, thus if:
% cv = [1 2 3 4]
% s1i = [5,14,13,11]
% s2i = [1,4,3,2]
% Then element 5 in a and element 1 in b are the closest to each other,
% with a distance of 1. 
% if numel(a)< numel(b)
%     temp = b;
%     b=a;
%     a=temp;
% end
  
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
    [cv,it] = min([abs(a-b(id));abs(b(iu)-a)]);
    ib = id+(it-1).*(iu-id);
% 
[cv,ix]=sort(cv);
s1i=ib(ix);
s2i=1:numel(a);
s2i=s2i(ix);
end
