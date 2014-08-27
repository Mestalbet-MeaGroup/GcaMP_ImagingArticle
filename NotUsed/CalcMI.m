function MutualInformation = CalcMI(traces,fr);
% Function which takes two timeseries and calculates the mutual information

if size(fr,2)>size(traces,2)
    subset  = randperm(size(fr,2),size(traces,2));
    fr=fr(:,subset);
end
if size(fr,2)<size(traces,2)
    subset  = randperm(size(traces,2),size(fr,2));
    traces=traces(:,subset);
end
fr=(fr-min(fr(:)))/(max(fr(:))-min(fr(:)));
traces=(traces-min(traces(:)))/(max(traces(:))-min(traces(:)));
% traces=otsu(traces,100);
% fr=otsu(fr,100);
indrow = double(uint16(fr(:))) + 1;
indcol = double(uint16(traces(:))) + 1; %// Should be the same size as indrow
jointHistogram = accumarray([indrow indcol], ones(length(indrow), 1));
jointProb = jointHistogram / length(indrow);
indNoZero = jointHistogram ~= 0;
jointProb1DNoZero = jointProb(indNoZero);
jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));

MutualInformation = entropy(fr)+entropy(traces)-jointEntropy;
end