function PlotBurstPropogationAnimation(e1,fr,num2hold)

% f = figure('Visible','off');
hold on;

%---Load Mea Image Data and Map---%
load('MeaMapPlot.mat','MeaImage','MeaMap');
MeaMap = rot90(MeaMap,-1);
overlay = otsu(imresize(MeaImage,1.0441),45);
overlay(overlay>1)=0;
overlay=bwareaopen(overlay, 500);
overlay = (overlay>=1);
%--Erase manually two spots---%
overlay(1750:1830,2705:2770)=0;
overlay(1600:1720,2745:2830)=0;
%-----------------------------%
centers = imfindcircles(overlay,[12,35],'ObjectPolarity','bright','method','TwoStage','Sensitivity',0.70);
[xpos,a]=sort(centers(:,2));
ypos=centers(a,1);

%---Plot Mea Image---%
% x = repmat([1:size(MeaImage,1)]',1,size(MeaImage,1));
% y = repmat(1:size(MeaImage,1),size(MeaImage,2),1);
% z = zeros(size(y));
% cdata = real2rgb(overlay,'gray',[1,0]);
% surface(x,y,z,cdata,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); axis('tight');axis('square');

%--Find Electrode Positions---%
mask = otsu(MeaImage,45);
mask(mask>1)=0;
centers = imfindcircles(mask,[12,20],'ObjectPolarity','bright','method','TwoStage','Sensitivity',0.90);
[xpos,a]=sort(centers(:,2));
ypos=centers(a,1);

%---Add Ground Electrodes to x and y pos vectors---%
xpos = [nan;xpos(1:14);nan;xpos(15:238);nan;xpos(239:end);nan];
ypos = [nan;ypos(1:14);nan;ypos(15:238);nan;ypos(239:end);nan];

%---Reshape xpos and ypos into 16x16 grid---%
xpos=reshape(xpos,16,16);
ypos=reshape(ypos,16,16);
% colors = jet(numel(e1));
colors = [158,255,0]./255;
[~,~,rank] = unique(fr);
hf=figure('color','white');
x=0;
y=0;
hold on
plot(centers(:,2),centers(:,1), 'o',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1]);
whitebg('black')
p = plot(x,y,'.','color',colors,'MarkerSize',25);
xlim([-100,3100]);
ylim([-100,3100]);

for i=1:numel(e1)
    [x1,y1] = find(MeaMap==e1(i));
    %     score=colors(rank(e1(i)),:);
    x = [x,xpos(x1,y1)];
    y = [y,ypos(x1,y1)];
    if i>num2hold
        x(i-num2hold)=nan;
        y(i-num2hold)=nan;
    end
    set(p,'XData',x)
    set(p,'YData',y)
    drawnow
    im = getframe;
    imwrite(im.cdata,'BurstAnimation.tif','WriteMode','append','compression','none');
end

end