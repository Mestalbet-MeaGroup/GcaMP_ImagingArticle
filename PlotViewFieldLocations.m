%---Viewfield Locations in RealSpace---%
load('MeaMapPlot.mat', 'MeaMap');
PosX = cellfun(@(x) (x.FramePosX), DataSet);
PosY = cellfun(@(x) (x.FramePosY), DataSet);
Layout = zeros(size(MeaMap,1)+1,size(MeaMap,2)+1);
for i=1:numel(PosX); Layout(PosX(i),PosY(i))=i; end
pcolor(Layout);
cmap = cbrewer('qual','Set1',9);
cmap = [[1 1 1];cmap];
colormap(gca,cmap);
colorbar;
set(gca,'BoxPlotAspectRatio',[1,1,1]);