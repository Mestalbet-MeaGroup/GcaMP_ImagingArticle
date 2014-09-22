function PlotImageAsSurface(mat);

xo = repmat(1:size(mat,2),size(mat,1),1);
yo = repmat([1:size(mat,1)]',1,size(mat,2));
zo = zeros(size(mat));
[rank,~,~] = unique(mat(:));
colors = [1,1,1;jet(numel(rank))];
cdatao = real2rgb(mat,colors,[0,max(rank)]);
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct');
colormap(colors);
colorbar;
end