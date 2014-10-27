function [prn,pra,prna,pran]=CalcPageRanks(MaxCosSim,ic);
a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
a2n=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);

% [~,~,~,a1,a_e2,~,ch,~]=ParseConnectionValues(w,1000);
% astros=[];
% for i=1:numel(unique(a_e2))
%     [test(i),F(i),~]=mode(a_e2);
%     astros = [astros,a1(ismember(a_e2,test(i)))];
%     a_e2(a_e2==test(i))=nan;
% end
% [~,~,~,a1,a_e2,~,ch,~]=ParseConnectionValuesWithVF(w,1000);
% 
% test(F<(mean(F)+std(F)))=[];
% astros(F<(mean(F)+std(F)))=[];
% 
% a2n = a2n(ismember(ic(1,:),a_e2),unique(a1));
% n2n = n2n(ismember(ic(1,:),a_e2),ismember(ic(1,:),a_e2));
% a2a = a2a(unique(a1),unique(a1));
% ic = ic(:,ismember(ic(1,:),a_e2));
% 
% MaxCosSim=zeros(size(n2n,1)+size(a2a,1),size(n2n,1)+size(a2a,1));
% MaxCosSim(1:size(n2n,1),1:size(n2n,2))=n2n+n2n'+eye(size(n2n));
% MaxCosSim(1:size(n2n,1),size(n2n,2)+1:end)=a2n;
% MaxCosSim(size(n2n,2)+1:end,1:size(n2n,1))=a2n';
% MaxCosSim(size(n2n,1)+1:end,size(n2n,1)+1:end)=a2a+a2a'+eye(size(a2a));
% 
%---Page Rank N2N---%
n2n = n2n+n2n'+eye(size(n2n));
n2n(logical(isnan(n2n)))=0;
prn = pagerank(n2n);
prn=prn';

%---Page Rank A2A---%
a2a = a2a+a2a'+eye(size(a2a));
a2a(logical(isnan(a2a)))=0;
pra = pagerank(a2a);
pra=pra';

%---Page Rank A2N---%
mat = zeros(size(n2n,1)+size(a2n,2),size(n2n,1)+size(a2n,2));
mat(1:size(ic,2),size(ic,2)+1:end)=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
mat=mat+mat'+eye(size(mat));
mat(isnan(mat))=0;

pr = pagerank(mat);
prna = pr(1:size(ic,2))';
pran = pr(size(ic,2)+1:end)';
end
