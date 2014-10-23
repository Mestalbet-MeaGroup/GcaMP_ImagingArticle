function [PeakTypeBurst,PeakTypeOther,on,nc,ClosestNear,ClosestFar]=ClassifyAstroPeaks(t,ic,traces,time,bs)
% Function which finds the peaks in the astro trace which are associated with bursts.

% [pi,~,~,ClosestValues,~]=FindPeaksNearestBursts(traces,time,bs,sort(t));
[~,~,ClosestValues,~,pi]=FindPeaksNearestBursts(traces,time,bs,sort(t));
%---Calculate the range for what to consider near to burst starts based on jPSTH---%
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

%---Assign the same time limits to each culture---%
sTime=1;
eTime=3;
%-------------------------------------------%

start =50;
after =50;
PeakTypeBurst=[];
PeakTypeOther=[];
on=[]; nc=[];
c=1; d=1;
for i=1:size(traces,2)
    if numel(ClosestValues{i})>1
        in  = find((ClosestValues{i}>0)&(ClosestValues{i}<=sTime)); %only look for traces which follow the burst start
        %         out = setdiff(1:numel(ClosestValues{i}),in);
        out = find(abs(ClosestValues{i})>eTime);
        OnIndex = pi{i}(in);
        ClosestNear{i} = ClosestValues{i}(in);
        ClosestFar{i} = ClosestValues{i}(out);
        notCloseOn = pi{i}(out);
        on = [time(OnIndex),on];
        nc = [time(notCloseOn),nc];
        
        
        for j=1:numel(OnIndex)
            offset = start-OnIndex(j);
            offsetE = OnIndex(j)+after-size(traces,1);
            if offsetE<0
                if offset<0
                    PeakTypeBurst(:,c) = traces(OnIndex(j)-start:OnIndex(j)+after,i);
                    c=c+1;
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    PeakTypeBurst(:,c) = temp(OnIndex(j)+offset-start+1:OnIndex(j)+after+offset+1);
                    c=c+1;
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                PeakTypeBurst(:,c) = temp(OnIndex(j)-start:OnIndex(j)+after);
                c=c+1;
            end
        end
        
        for j=1:numel(notCloseOn)
            offset = start-notCloseOn(j);
            offsetE = notCloseOn(j)+after-size(traces,1);
            if offsetE<0
                if offset<0
                    PeakTypeOther(:,d) = traces(notCloseOn(j)-start:notCloseOn(j)+after,i);
                    d=d+1;
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    PeakTypeOther(:,d) = temp(notCloseOn(j)+offset-start+1:notCloseOn(j)+after+offset+1);
                    d=d+1;
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                PeakTypeOther(:,d) = temp(notCloseOn(j)-start:notCloseOn(j)+after);
                d=d+1;
            end
        end
        
    end
end

end