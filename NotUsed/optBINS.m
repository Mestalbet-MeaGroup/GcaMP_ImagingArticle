function optM = optBINS(data,minM,maxM)
% Taken from:
% @ARTICLE{2006physics...5197K,
%    author = {{Knuth}, K.~H.},
%     title = "{Optimal Data-Based Binning for Histograms}",
%   journal = {ArXiv Physics e-prints},
%    eprint = {physics/0605197},
%  keywords = {Physics - Data Analysis, Statistics and Probability, Mathematics - Probability, Mathematics - Statistics Theory, Physics - Computational Physics},
%      year = 2006,
%     month = may,
%    adsurl = {http://adsabs.harvard.edu/abs/2006physics...5197K},
%   adsnote = {Provided by the SAO/NASA Astrophysics Data System}

if size(data)>2 | size(data,1)>1
    error('data dimensions must be (1,N)');
end
N = size(data,2);
logp = zeros(1,maxM);
% Simply loop through the different numbers of bins
% and compute the posterior probability for each.
for M = minM:maxM
    n = hist(data,M); % Bin the data (equal width bins here)
    logp(M) = N*log(M) + gammaln(M/2) - gammaln(N+M/2) - M*gammaln(1/2) + sum(gammaln(n+0.5));
end
[~, optM] = max(logp(minM,maxM));
end