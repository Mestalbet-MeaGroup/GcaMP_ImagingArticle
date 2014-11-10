function [neuro,astro,n2a,a2n] = CalcLinkToNetwork(DataSet,k)

mtr= nanmean(DataSet{k}.dfTraces,2);
gfr= DataSet{k}.GFR;
[num,~] = max([size(DataSet{k}.dfTraces,2),size(DataSet{k}.FR,2)]);
for i=1:num
    if i<=size(DataSet{k}.dfTraces,2)
        a2n(i) = CalcMaxCorr(DataSet{k}.dfTraces(:,i),gfr);
        astro(i) = CalcMaxCorr(DataSet{k}.dfTraces(:,i),mtr);
    end
    if i<=size(DataSet{k}.FR,2)
        neuro(i) = CalcMaxCorr(DataSet{k}.FR(:,i),gfr);
        n2a(i) = CalcMaxCorr(DataSet{k}.FR(:,i),mtr);
    end
end
