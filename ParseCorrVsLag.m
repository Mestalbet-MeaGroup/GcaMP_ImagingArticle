function vals = ParseCorrVsLag(which);
CS  = matfile('F:\CosSimTemp.mat');
csi = CS.CosSim(1,which);
csi=csi{1};
% csi = reshape(csi,333*333,10000);
c=1;
% for i=1:size(csi,1)
%     for j=1:size(csi,2)
%         vals(c,1) = csi(i,j);
%         vals(c,2) = j;
%         c=c+1;
%     end
% end
% end

%---Perhaps a faster way---%
CS  = matfile('F:\CosSimTemp.mat');
csi = CS.CosSim(1,which);
csi=csi{1};
csj = sparse([],[],[],size(csi,1)*size(csi,2),size(csi,3));
c=1;
for i=1:size(csi,1)-1
    for j=i+1:size(csi,2)
        csj(c,:)=squeeze(csi(i,j,:));
        c=c+1;
    end
end

clear csi; 
[xi,yi]=ind2sub(size(csi),1:numel(csi));
[xj,yj]=ind2sub(size(csj),1:numel(csi));
