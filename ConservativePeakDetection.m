function p = ConservativePeakDetection(traces)
around = 100;
fun = @(block_struct) max(block_struct.data);
% for i=1:size(traces,2)
    [~,idx,~] = deleteoutliers(traces,0.01,0);
    test = zeros(max(idx),1);
    test(idx)=traces(idx);
    res = blockproc(test,[150,1],fun);
    res(res==0)=[];
    bs = ismember(traces,res);
    bs = find(bs);
    tr = padarray(traces,[around,0],nan,'both');
    bs = arrayfun(@(x) find(tr==nanmax(tr(x-around:x+around)),1,'First'),bs+around)-around;
    p=unique(bs);
%     p{i}=bs;
% end