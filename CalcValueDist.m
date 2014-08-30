function [distance,score] = CalcValueDist(mat,ic,mask,pos,options);
% Calculates the distance in real space between elements in mat
% Options; 'a2n', 'n2n' or 'a2a'
%Usage:
% N2N:  [distance,score] = CalcValueDist(n2n,ic,[],[],'n2n');
% A2A:  [distance,score] = CalcValueDist(a2a,[],mask,[],'a2a');
% A2N:  [distance,score] = CalcValueDist(a2n,ic,mask,173,'a2n');
n = size(mat,1);
m = size(mat,2);

switch options
    case 'a2a'
        if m~=n
            error('You chose A2A but the matrix is not symmetric, not the same number of Astros');
        end
        distmat = regiondist(~mask);
        distmat(distmat==0)=nan;
        distmat = triu(distmat,1);
        distmat(isnan(distmat))=-1;
        distmat(distmat==0)=nan;
        distmat(distmat==-1)=0;
        mat(mat==0)=nan;
        mat = triu(mat,1);
        mat(isnan(mat))=-1;
        mat(mat==0)=nan;
        mat(mat==-1)=0;
        distance = distmat(:);
        distance(isnan(distance))=[];
        score = mat(:);
        score(isnan(score))=[];
    case 'a2n'
        MeaImage = [];
        MeaMap=[];
        load('MeaMapPlot.mat','MeaImage','MeaMap');
        %         MeaMap = rot90(MeaMap,-1);
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
        
        xpos=round(xpos);
        ypos=round(ypos);
        xpos = [nan;xpos(1:14);nan;xpos(15:238);nan;xpos(239:end);nan];
        ypos = [nan;ypos(1:14);nan;ypos(15:238);nan;ypos(239:end);nan];
        xpos=reshape(xpos,16,16);
        ypos=reshape(ypos,16,16);
        
        %         mask=imresize(mask,0.6472); % There is a problem with imresize - it deletes thin bridges.
        %---Resize image without loosing small bridges---%
        %         mask = imwarp(double(mask)+100,affine2d([0.6472, 0, 0;0, 0.6472,0; 0,0,1]),'Interp','linear');
        %         mask(mask==101)=1;
        %         mask(mask<100)=1;
        %         mask(mask~=1)=0;
        % %-------------------------------------------------%
        [sitex,sitey]=find(MeaMap==pos);
        centframe = [round(size(mask,1)/2),round(size(mask,2)/2)];
        x0 = xpos(sitex,sitey)-centframe(1);
        y0 = ypos(sitex,sitey)-centframe(2); %convert center of frame location to corner
      
%         ROI = regionprops(~mask,'PixelList');
        if m==size(ic,2)
            numAstros=n;
        else
            numAstros=m;
        end
        
%         if numel(ROI)~=numAstros
%             error('Wrong number of ROIs detected.');
%         end
        distmat=nan(size(mat));
        labeled = bwlabel(~mask);
        if max(labeled(:))~=numAstros
            error('Wrong number of ROIs detected.');
        end
        T = maketform('affine',[0.6472, 0, 0;0, 0.6472,0; 0,0,1]); %The transform to match scales of mea layout and mask. 
        for i=1:max(labeled(:));
            [u,v] = find(labeled==i);
            [x,y] = tformfwd(T,u,v);
            temp=CalcDistanceROIElectrodeCenter(MeaMap,ic,xpos,ypos,overlay,false(size(mask)),[x,y],x0,y0);
            distmat(:,i)=temp;
        end
        distance=distmat(:);
        score = mat(:);
    case 'n2n'
        if m~=n
            error('You chose N2N but the matrix is not symmetric, not the same number of neuros');
        end
        MeaImage=[];
        MeaMap=[];
        load('MeaMapPlot.mat','MeaImage','MeaMap');
        %         MeaMap = rot90(MeaMap,-1);  %check to make sure MeaMap doesn't need to rotated
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
        
        combs=VChooseK(1:n,2);
        e1=ic(1,combs(:,1));
        e2=ic(1,combs(:,2));
        distance=zeros(numel(e1));
        score=zeros(numel(e1));
        for i=1:numel(e1)
            [x1,y1] = find(MeaMap==e1(i));
            [x2,y2] = find(MeaMap==e2(i));
            distance(i) = sqrt((xpos(x1,y1)-xpos(x2,y2))^2+ (ypos(x1,y1)-ypos(x2,y2))^2);
            score(i)=mat(combs(i,1),combs(i,2));
        end
end

end

