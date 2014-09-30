function [radius, Corr]=CalcCorrVsRadius(mask,corrmat,ic,ch);

%---Get VF Electrodes---%
load('MeaMapPlot.mat','MeaMap');
[~,~,channels] = FindElecsinVF(ch,MeaMap);
channels = reshape(channels,3,3)';
[~,col]=find(channels==ch);
channels=fliplr(channels(col,:));

%---Find which ROIs intersect with circle---%
pixels = regionprops(~mask,'PixelIdxList');
Corr = zeros(size(corrmat,2),150);
load('Cult6_ElecLocsViewField.mat');
circs = regionprops(~elecs,'Centroid','PixelIdxList');
circs = circs(arrayfun(@(x) numel(x.PixelIdxList)>100,circs));
parfor r=1:150
    Corr(:,r) = CalcCvsR(corrmat,ic,mask,pixels,circs,channels,ch,r);
end
radius = 1:150;
% Corr = nanmean(Corr,1);

    
end