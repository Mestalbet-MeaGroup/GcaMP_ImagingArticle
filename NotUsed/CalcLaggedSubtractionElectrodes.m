vec = ConvertIC2Samora(ic);
list=unique(vec);
temp = nchoosek(1:numel(list),2);
set1=list(temp(:,1));
set2=list(temp(:,2));

for i=1:numel(list)
    fr(i,:) = histc(t(vec==list(i)),0:100:max(t));
end


mat=zeros(numel(list),numel(list));
maxlag=100;
for j=1:maxlag
    for i =1:numel(set1)
        mat(temp(i,1),temp(i,2),j) = 1 -(numel([setdiff(t(vec==set1(i))+j,t(vec==set2(i))),setdiff(t(vec==set2(i)),t(vec==set1(i))+j)])./(sum(vec==set1(i))+sum(vec==set2(i))));
    end
end

