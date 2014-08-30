function distmat=CalcDistanceROIElectrodeCenter(MeaMap,ic,xpos,ypos,overlay,tempmask,pixs,x0,y0)
numc = size(ic,2);
distmat = nan(numc,1);
for ii=1:size(pixs,1);
    tempmask(round(pixs(ii,1)),round(pixs(ii,2)))=1;
end
for j=1:numc
    [sitex,sitey]=find(MeaMap==ic(1,j));
    elecx = xpos(sitex,sitey);
    elecy = ypos(sitex,sitey);
    tempoverlay=false(size(overlay));
    tempoverlay(elecx,elecy)=1;
    
    %     tempmask=(~mask);
    %     tempoverlay=false(size(overlay));
    %     for j=1:numc
    %         [sitex,sitey]=find(MeaMap==ic(1,j));
    %         elecx = xpos(sitex,sitey);
    %         elecy = ypos(sitex,sitey);
    %         tempoverlay(elecx,elecy)=1;
    %     end
    
    %---For cases when viewfield goes beyond where there are electrodes.Centered on border electrode-%
    if ((x0+size(tempmask,1)-1)>size(tempoverlay,1))
        tempoverlay= padarray(tempoverlay,[(x0+size(tempmask,1)-1)-size(tempoverlay,1),0],0,'post');
    end
    
    if ((y0+size(tempmask,2)-1)>size(tempoverlay,2))
        tempoverlay = padarray(tempoverlay,[0,(y0+size(tempmask,2)-1)-size(tempoverlay,2)],0,'post');
    end
    
    if (x0<0)
        tempoverlay= padarray(tempoverlay,[abs(x0),0],0,'pre');
        x0 = x0+abs(x0)+1;
    end
    
    if (y0<0)
        tempoverlay= padarray(tempoverlay,[0,abs(y0)],0,'pre');
        y0 = y0+abs(y0)+1;
    end
    %-----------------------------------------------------------------------------------------------%
    
    tempoverlay(x0:x0+size(tempmask,1)-1,y0:y0+size(tempmask,2)-1)=...
        tempoverlay(x0:x0+size(tempmask,1)-1,y0:y0+size(tempmask,2)-1)+(tempmask);
    %     temp = regiondist(tempoverlay);
    %     if size(temp,1)>2
    %         printf('Problem');
    %     end
    %     if size(temp,1)<2
    %         distmat(j,1)=0;
    %     else
    %         distmat(j,1)=temp(1,2);
    %     end
    label = bwlabel(tempoverlay);
    if sum(label(:)==2)>0
        dt = bwdist(label == 1);
        [distmat(j,1),~] = min(dt(label==2));
    else
        distmat(j,1)=0;
    end
    
end