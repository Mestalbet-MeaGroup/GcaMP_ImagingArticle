function [PeakTypeBurst,PeakTypeOther]=GroupTracesNearBurstStart(t,ic,traces,time,bs)
% Function which finds the peaks in the astro trace which are associated with bursts.


%---Calculate the Range for near burst satrt peaks---%
% [PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);
% clusters = kmeans(PairWisePSTH',5);
% [~,ix]=max([sum(clusters==1),sum(clusters==2),sum(clusters==3),sum(clusters==4),sum(clusters==5)]);
% clusters(clusters~=ix)=0;
% clusters(clusters==ix)=1;
% [sIDX,eIDX]=initfin(~clusters');
% for i=1:numel(sIDX)
%     sc(i) = nansum(nansum(PairWisePSTH(:,sIDX(1):eIDX(1)),2));
% end
% sc(PairWiseLags(sIDX-1)<0)=nan;
% [~,w]=nanmax(sc);
% we = find(diff(eIDX(w:end))<5,1,'last')+w;
% sTime = PairWiseLags(sIDX(w)-1)/1000;
% if ~isempty(we)
%     eTime = PairWiseLags(eIDX(we))/1000;
% else
%     eTime = PairWiseLags(eIDX(w))/1000;
% end
% sTime=0.7;
eTime=1.5;
%-------------------------------------------%

start = 1;
after = 49;
PeakTypeBurst=[];
PeakTypeOther=[];
c=1;
for i=1:size(traces,2)
    for j=1:numel(bs)
        OnIndex = find(((time>=(bs(j)/12000))&(time<(bs(j)/12000+eTime))),1,'First');
        if ~isempty(OnIndex)
            offset = start-OnIndex;
            offsetE = OnIndex+after-size(traces,1);
            if offsetE<0
                if offset<0
                    temp = traces(OnIndex-start:OnIndex+after,i);
                    [~,idx,~] = deleteoutliers(temp,0.01);
                    if numel(idx)>1
                        PeakTypeBurst(:,c) = temp;
                        c=c+1;
                    end
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    [~,idx,~] = deleteoutliers(temp,0.01);
                    if numel(idx)>1
                        PeakTypeBurst(:,c) = temp(OnIndex+offset-start+1:OnIndex+after+offset+1);
                        c=c+1;
                    end
                end
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                [~,idx,~] = deleteoutliers(temp,0.01);
                if numel(idx)>1
                    PeakTypeBurst(:,c) = temp(OnIndex-start:OnIndex+after);
                    c=c+1;
                end
            end
        end
    end
end
c=1;
for i=1:size(traces,2)
    for j=1:numel(bs)
        notCloseOn = find(time>=(bs(j)/12000+eTime),1,'First');
        if ~isempty(notCloseOn)
            offset = start-notCloseOn;
            offsetE = notCloseOn+after-size(traces,1);
            if offsetE<0
                if offset<0
                    temp = traces(notCloseOn-start:notCloseOn+after,i);
                    [~,idx,~] = deleteoutliers(temp,0.01);
                    if numel(idx)>1
                        PeakTypeOther(:,c) = temp;
                        c=c+1;
                    end
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    [~,idx,~] = deleteoutliers(temp,0.01);
                    if numel(idx)>1
                        PeakTypeOther(:,c) = temp(notCloseOn+offset-start+1:notCloseOn+after+offset+1);
                        c=c+1;
                    end
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                [~,idx,~] = deleteoutliers(temp,0.01);
                if numel(idx)>1
                    PeakTypeOther(:,c) = temp(notCloseOn-start:notCloseOn+after);
                    c=c+1;
                end
            end
        end
    end
end
end