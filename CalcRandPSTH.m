function psth = CalcRandPSTH(t,ic,starts,ends,ttime,upto,numlag)
%calculates psth for randomly shifted interval spikes
permt=zeros(size(t));
pre= 5000; %ms
post = pre;
bin = 100; %ms
psth = zeros(upto,numlag);

for i=1:size(ic,2)
    temp = ShuffelIntervals([min(t(ic(3,i):ic(4,i))),sort(t(ic(3,i):ic(4,i)))]);
    permt(ic(3,i):ic(4,i))=temp(2:end);
%     temp=permt(ic(3,i));
%     permt(ic(3,i))=permt(ic(4,i));
%     permt(ic(4,i))=temp;
end
[permt,ic]=CutSortChannel2(permt,ic,floor(min(cellfun(@(x)min(x),ttime))-5).*12000,max(permt));
fun = @(s,e) permt(s:e);
tt = arrayfun(fun, starts,ends,'UniformOutput',false);
parfor c=1:upto
    ps = mpsth(tt{c}./12000,ttime{c},'tb',0,'fr',1,'pre',pre,'post',post,'binsz',bin)./numel(ttime{c});
    psth(c,:) = ps;
end
end