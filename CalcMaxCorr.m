function Corr=CalcMaxCorr(a,b)
a=zscore(a);
b=zscore(b);
len = length(a);
auto_a = xcorr(a,'unbiased');
auto_b = xcorr(b,'unbiased');
[Corr,~] = xcorr(a,b,1000,'unbiased');
Corr = Corr ./ sqrt(auto_a(len) .* auto_b(len)); % normalize by values at zero lag
Corr  = max(Corr);
end