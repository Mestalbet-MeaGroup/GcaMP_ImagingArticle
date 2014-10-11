function [u,v] = ConvertResizedCoordinates(mask,overlay,x,y)
% Function which takes x,y and the ROI mask with its resized mask
% (overlay), and calculates the transformed coordinates.
rowsinA = size(mask,2);
rowsinB = size(overlay,2);
r = rowsinB / rowsinA;

colsinA = size(mask,1);
colsinB = size(overlay,1);
q = colsinB / colsinA;

u = r * (x - 0.5) + 0.5;
v = q * (y - 0.5) + 0.5;
end