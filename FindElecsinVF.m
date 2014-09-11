function [xind,yind] = FindElecsinVF(ch,MeaMap)
% Calculates the index of the electrodes surrounding the electrode in the
% view field center. 
[i,j]=find(MeaMap==ch);
[Iadj , Radj, Nfound ] = neighbourND(sub2ind(size(MeaMap),i,j), size(MeaMap));
[xind,yind] = ind2sub(size(MeaMap),Iadj);
end
