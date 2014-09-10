function [counts,bins]=CalcHist(x)
    [counts,~,bins] = histomex(x',nsOPTBINS(x')); 
    bins=bins{1};
end