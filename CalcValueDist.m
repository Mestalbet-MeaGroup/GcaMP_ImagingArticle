function [distance,score] = CalcValueDist(mat,ic,mask,options);
% Calculates the distance in real space between elements in mat
% Options; 'a2n', 'n2n' or 'a2a'
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
        
    case 'n2n'
        load('MeaMapPlot.mat','MeaImage','MeaMap');
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
        for i=1:numel(e1)
            [x1,y1] = find(MeaMap==e1(i));
            [x2,y2] = find(MeaMap==e2(i));
            distance(i) = sqrt((xpos(x1,y1)-xpos(x2,y2))^2+ (ypos(x1,y1)-ypos(x2,y2))^2);
            score(i)=mat(combs(i,1),combs(i,2));
        end       
end

end
