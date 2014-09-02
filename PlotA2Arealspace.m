function PlotA2Arealspace(mask,connmat,number);
% Function which accepts the mask, connection matrix and the number of
% lines to plot and produces a plot of the connections between astrocytes
% in real space.

connmat1= triu(connmat,1);
for i=1:number
    [avals(i),loc(i)]=nanmax(connmat1(:)); 
    connmat1(loc(i))=nan; 
end
[a1,a2]=ind2sub(size(connmat1),loc);
% [centrality,index]=CalcEigenCentrality(connmat,min(avals));
[centrality,index]=CalcPageRank(connmat);
centrality = round(abs(centrality).*10000)./10000;
[~,~,noderank] = unique(centrality,'stable');
nodecolors = [0,0,0;flip(cbrewer('seq','Purples',numel(unique(noderank)-2)),1)];
pixels = regionprops(~mask,'PixelIdxList');

image = zeros(size(mask));
for i=1:numel(centrality)
    image(pixels(index(i)).PixelIdxList)=noderank(i);
end

cdatao = real2rgb(image,nodecolors,[0,max(noderank)]);
xo = repmat(1:size(mask,2),size(mask,1),1);
yo = repmat([1:size(mask,1)]',1,size(mask,2));
zo = zeros(size(mask));
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); axis('tight');

% cvec = [0,0,0;nodecolors(2:end,:)];
colormap([0,0,0;flip(nodecolors(2:end,:),1)]);
c = colorbar('Location','northoutside');
cticks = linspace(0,1,size(nodecolors,1));
set(c,'Ticks',[0,cticks(2),cticks(end-1)],'TickLabels',[0,1,max(noderank)]);
ylabel(c,'Page Rank')
caxis(caxis);


centers  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~mask,'Centroid'),'uniformoutput',0));
[avals,rank]=sort(round(avals*10000)/10000);
a1=a1(rank);
a2=a2(rank);
[~,~,avals_index] = unique(avals);
rank = rank(avals_index);
xposa1 = (centers(a1,1))';
yposa1 = (centers(a1,2))';
xposa2 = (centers(a2,1))';
yposa2 = (centers(a2,2))';
colors = jet(max(rank)-min(rank)+1);
colors = colors(rank-min(rank)+1, :);
%--Transparent Lines---%
for i=1:numel(a1)
    l = patchline([xposa1(i),xposa2(i)],[yposa1(i),yposa2(i)],'linestyle','-','edgecolor',colors(i,:),'linewidth',1,'edgealpha',0.5);
end

%--Solid Lines---%
% x=[xposa1;xposa2];
% xt(:,1:2:size(x,2)*2)=x;
% xt(:,xt(1,:)==0)=nan;
% y=[yposa1;yposa2];
% yt(:,1:2:size(y,2)*2)=y;
% yt(:,yt(1,:)==0)=nan;
% counter=1;
% for i=1:2:size(xt,2)-2; 
%     vecx(counter,:)=[xt(1,i),xt(2,i),xt(1,i+1),xt(1,i+2),xt(2,i+2)]; 
%     vecy(counter,:)=[yt(1,i),yt(2,i),yt(1,i+1),yt(1,i+2),yt(2,i+2)]; 
%     counter=counter+1; 
% end
% l = line(vecx',vecy','linewidth',1);
% hline = findobj(gcf, 'type', 'line');
% for i=1:size(hline)
%     set(hline(i),'color',colors(i,:));
% end
%--------------%

a1labels =cellfun(@(x) num2str(x), num2cell(a1),'UniformOutput',false);
text(xposa1,yposa1,a1labels,'Clipping', 'on','hittest','off','color','r');

ha  = axes('visible','off');
colormap(ha,jet(max(rank)-min(rank)+1));
c1 = colorbar('Location','eastoutside');
ctick = 0:(1/numel(avals))/2:1;
ctick=ctick(2:2:end-1);
set(c1,'Ticks',[ctick(1),ctick(end)],'TickLabels',[avals(1),avals(end)]);
ylabel(c1,'Correlation')
end