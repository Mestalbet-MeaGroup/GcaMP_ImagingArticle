function corr = CalcCvsR(corrmat,ic,mask,pixels,circs,channels,ch,r)
% Function which calculates correlation as a function of ROI within radius
% of electrode. TO be used to with CalcCorrVsRadius.m
corr = zeros(size(corrmat,2),1);
for i=1:size(corrmat,2)
    tempmask = false(size(mask));
    tempmask(pixels(i).PixelIdxList)=true;
    %---Create Circles of Radius R around electrodes---%
    cc=1:size(mask,2);
    rr=[1:size(mask,1)]';
    c=0;
    for j=1:3
        %---Select Electrode---%
        cx1=circs(j).Centroid(2); 
        cy1=circs(j).Centroid(1);
        f=@(xx,yy) (xx-cx1).^2+(yy-cy1).^2 <=r^2; % | (xx-cx2).^2+(yy-cy2).^2 <=r^2 | (xx-cx3).^2+(yy-cy3).^2 <=r^2;
        %---Place circle  at Electrode Center---%
        tempcirc=logical(bsxfun(f,rr,cc));
        %---Check if there overlap---%
        overlap = tempmask&tempcirc;
        if sum(overlap(:))>0;
            c=c+1;
            elec = channels(find(channels==ch)+(j-2));
            corr(i)= corr(i)+corrmat(ic(1,:)==elec,i);
        end
    end
    corr(i)=(corr(i)/c)/nanmean(corrmat(:,i));
end
end