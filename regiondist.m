function d = regiondist(img)
label = bwlabel(img);
n = max(label(:));
d = zeros(n);
for ii = 1:n-1
  dt = bwdist(label == ii);
  for jj = ii+1:n
    reg = (label == jj);
    [mindist,~] = min(dt(reg));
    d(ii,jj) = mindist;
  end
end
d=d+d';
end