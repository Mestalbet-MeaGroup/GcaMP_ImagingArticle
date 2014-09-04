function CosSim = CalcCosineSimilarity(fr,traces,maxlag)
% Function which returns the lagged cosine similarity (where input traces have
% means subtracted) equivelant to correlation coefficient.

% CosSim = zeros(size(fr,2)+size(traces,2),size(fr,2)+size(traces,2),maxlag);
lim=floor(maxlag/2)+1;
CosSimNeg=zeros(size(fr,2)+size(traces,2),size(fr,2)+size(traces,2),lim);
% CosSimPos=CosSimNeg;

parfor j=1:lim
        ifr = (lim+1-j):size(fr,1);
        itrace = 1:(size(traces,1)-lim+j);
        r = [fr(ifr,:),traces(itrace,:)]';
        r=r-repmat(mean(r,2),1,size(r,2));
        norm_r = sqrt(sum(abs(r).^2,2));
        r= r./repmat(norm_r,1,size(r,2));
        tr = transpose(r);    
        CosSimNeg(:,:,j) = (r*tr);
       
        r = [fr(itrace,:),traces(ifr,:)]';
        r=r-repmat(mean(r,2),1,size(r,2));
        norm_r = sqrt(sum(abs(r).^2,2));
        r= r./repmat(norm_r,1,size(r,2));
        tr = transpose(r);    
        CosSimPos(:,:,j) = (r*tr);
end
% CosSim=cat(3,CosSimNeg,flip(CosSimNeg(:,:,1:end-1),3));
% CosSimPos(:,:,maxlag/2+1:end)=[];
% CosSimNeg(:,:,1:maxlag/2)=[];
% CosSimNeg = flip(CosSimNeg,3);
CosSim=cat(3,CosSimNeg,flip(CosSimPos(:,:,1:end-1),3));
end