function PlotHistOverExpected(data,optM)
[temphist,bins]=hist3(data,optM);
counts = repmat(sum(temphist,2),[1,optM(1)]);
counts = counts/sum(counts(:,1));
counts2 = repmat(sum(temphist,1),[optM(2),1]);
counts2 = counts2/sum(counts2(1,:));
expected = counts.*counts2;
expected = expected./trapz(bins{2},trapz(bins{1},expected));

[temphist,bins]=hist3(data,optM);
normhist = temphist/trapz(bins{2},trapz(bins{1},temphist));

norm2expected = normhist-expected;
if mod(numel(unique(norm2expected(:))),2)>0
    colors = cbrewer('div','RdBu',numel(unique(norm2expected(:)))+1);
    colors = [colors(1:floor(end/2),:);[1,1,1];colors(floor(end/2)+1:end,:)];
else
    colors = cbrewer('div','RdBu',numel(unique(norm2expected(:))));
    colors = [colors(1:floor(end/2),:);[1,1,1];colors(floor(end/2)+1:end,:)];
end
colors = colors(end:-1:1,:);
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);

surface(xo,yo,round(norm2expected.*10000)/10000);
colormap(colors);
caxis([-max(abs(norm2expected(:))),max(abs(norm2expected(:)))]);
hCbar = colorbar;
axis tight;
end