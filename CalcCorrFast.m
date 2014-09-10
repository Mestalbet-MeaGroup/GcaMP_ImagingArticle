function c1=CalcCorrFast(x,maxlag)

combs = VChooseKR(1:size(x,1),2);
e1 = combs(:,1);
e2 = combs(:,2);
c=zeros(size(x,1)*size(x,1),maxlag*2+1);
linidx= sub2ind([size(x,1),size(x,1)],e1,e2);
%---Parralel Option 1---%
% spmd
%     temp=zeros(numel(e1),maxlag*2+1,codistributor1d(1));
%     for i=drange(1:numel(e1))
%         if abs(e1(i)-e2(i))>0
%             temp(i,:)=xcorr(x(e1(i),:),x(e2(i),:),maxlag,'coeff');
%         end
%     end
% end
% temp = gather(temp);
% c(linidx,:)=temp;

%---Parralel Option 2---%
temp=zeros(numel(e1),maxlag*2+1);
parfor i=1:numel(e1)
    if abs(e1(i)-e2(i))>0
        temp(i,:)=xcorr(x(e1(i),:),x(e2(i),:),maxlag,'coeff');
    end
end
c(linidx,:)=temp;

%---Serial Option---%
% for i=1:numel(e1)
%     if abs(e1(i)-e2(i))>0
%         c(linidx(i),:)=xcorr(x(e1(i),:),x(e2(i),:),maxlag,'coeff');
%     end
% end
%------%
c1 = reshape(c,[size(x,1),size(x,1),maxlag*2+1]);
end