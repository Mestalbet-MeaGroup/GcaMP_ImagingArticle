function [h,e,c] = PlotHistWithErrors(data,M,thickness,color,edges,centers)
% M is the number of bins
% if edges is empty, calculates the histogram bin centers and counts per
% bin. Otherwise, uses supplied edges to calculate bin counts.

e=[]; c=centers;
[D, N] = size(data);
if isempty(edges)
    [n, binwidth,centers,e] = histogram(data,M);
    c = centers{1}; % determine the bin counts
else
    binwidth = mean(diff(edges(2:end-1)));
    factors = ones(1,D);
    counts = zeros(1,M);
    for d = 1:D-1
        factors(d+1:D) = factors(d+1:D)*M(d);
    end
    
    for n = 1:N
        bin = zeros(1,D);
        for d = 1:D
            for m = 1:M(d)  %loop through all the bins in this dimension
                if (data(d,n) <= edges(d,m))
                    bin(d) = m;
                    break;
                end
            end
        end
        index = factors * (bin-1)' + 1;        % now add a count to that bin
        if index~=0
            counts(index) = counts(index)+1;
        end
    end
    n=counts;
end
p = ((n+0.5)/(N+M/2)/binwidth);                                % determine the mean bin probabilities
s = (sqrt(p .* (N-n+(M-1)/2)/((N+M/2+1)*(N+M/2)))/binwidth);   % determine the std devs
h = bar(c, p,thickness,'FaceColor',color);
for m = 1:M
    line([c(m); c(m)], [p(m); p(m)-s(m)], 'Color', 'w');
    line([c(m)- binwidth/6; c(m)+ binwidth/6], [p(m)-s(m); p(m)-s(m)], 'Color', 'w');
    line([c(m); c(m)], [p(m); p(m)+s(m)], 'Color', 'k');
    line([c(m)- binwidth/6; c(m)+ binwidth/6], [p(m)+s(m); p(m)+s(m)], 'Color', 'k');
end
end