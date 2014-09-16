function [PeakTypeBurst,PeakTypeOther,trCutouts]=ClassifyAstroPeaks(t,ic,traces,time,bs,be)
% Function which finds the peaks in the astro trace which are associated with bursts.

[pi,oi,ofi,ClosestValues]=FindPeaksNearestBursts(traces,time,bs);

% OffIndex  = ofi{j}(clusts==ix);
% peakIndex = pi{j}(clusts==ix);
PeakTypeBurst=[];
PeakTypeOther=[];
subplot(1,2,1)
hold on;
for i=1:size(traces,2)
    if numel(ClosestValues{i})>5
        clusts = kmeans(ClosestValues{i}',3);
        [~,ix] = min([nanmean(ClosestValues{i}(clusts==1)),nanmean(ClosestValues{i}(clusts==2)),nanmean(ClosestValues{i}(clusts==3))]);
        OnIndex = oi{i}(clusts==ix);
        for j=1:numel(OnIndex)
            offset = 10-OnIndex(j);
            offsetE = OnIndex(j)+60-size(traces,1);
            if offsetE<0
                if offset<0
                    plot(traces(OnIndex(j)-10:OnIndex(j)+60,i));
                    PeakTypeBurst(:,j) = traces(OnIndex(j)-10:OnIndex(j)+60,i);
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    plot(temp(OnIndex(j)+offset-10+1:OnIndex(j)+60+offset+1));
                    PeakTypeBurst(:,j) = temp(OnIndex(j)+offset-10+1:OnIndex(j)+60+offset+1);
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                plot(temp(OnIndex(j)-10:OnIndex(j)+60));
                PeakTypeBurst(:,j) = temp(OnIndex(j)-10:OnIndex(j)+60);
            end
        end
    end
end


subplot(1,2,2)
hold on;
for i=1:size(traces,2)
    if numel(ClosestValues{i})>5
        clusts = kmeans(ClosestValues{i}',3);
        [~,ix] = min([nanmean(ClosestValues{i}(clusts==1)),nanmean(ClosestValues{i}(clusts==2)),nanmean(ClosestValues{i}(clusts==3))]);
        notCloseOn = oi{i}(clusts~=ix);
        for j=1:numel(notCloseOn)
            offset = 10-notCloseOn(j);
            offsetE = notCloseOn(j)+60-size(traces,1);
            if offsetE<0
                if offset<0
                    plot(traces(notCloseOn(j)-10:notCloseOn(j)+60,i));
                    PeakTypeOther(:,j) = traces(notCloseOn(j)-10:notCloseOn(j)+60,i);
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    plot(temp(notCloseOn(j)+offset-10+1:notCloseOn(j)+60+offset+1));
                    PeakTypeOther(:,j) = temp(notCloseOn(j)+offset-10+1:notCloseOn(j)+60+offset+1);
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                plot(temp(notCloseOn(j)-10:notCloseOn(j)+60));
                PeakTypeOther(:,j) = temp(notCloseOn(j)-10:notCloseOn(j)+60);
            end
        end
    end
end

%---Method 2---%
[PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);
clusters = kmeans(PairWisePSTH',2);
[~,ix]=max([sum(clusters==1),sum(clusters==2)]);
clusters(clusters~=ix)=0;
clusters(clusters==ix)=1;

[sIDX,eIDX]=initfin(~clusters');
sTime = PairWiseLags(sIDX(1))/1000;
eTime = PairWiseLags(eIDX(1))/1000;
bs = bs./12000;
c=1; 
trCutouts=nan(200,size(traces,2),numel(bs));
for i=1:numel(bs)
    if (bs(i)+sTime)>=min(time)
        if (bs(i)+eTime)<=max(time)
            temp = traces( (time>(bs(i)-sTime))&(time<(bs(i)+eTime*3)),:);
            trCutouts(1:size(temp,1),1:size(temp,2),c) = temp;
            c=c+1;
        end
    end
end
trCutouts = reshape(trCutouts,[200,size(trCutouts,2)*size(trCutouts,3)]);
end