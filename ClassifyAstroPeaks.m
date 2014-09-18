function [PeakTypeBurst,PeakTypeOther,on,nc]=ClassifyAstroPeaks(t,ic,traces,time,bs)
% Function which finds the peaks in the astro trace which are associated with bursts.

[pi,~,~,ClosestValues]=FindPeaksNearestBursts(traces,time,bs);
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
sTime=0.7;
% eTime=1.5;
%-------------------------------------------%

start =30;
after =30;
PeakTypeBurst=[];
PeakTypeOther=[];
on=[]; nc=[];

% figure('visible','off','Position', [0 0 screen_size(3) screen_size(4) ])
% subplot(1,2,1)
% hold on;
for i=1:size(traces,2)
    if numel(ClosestValues{i})>1
        in  = find(ClosestValues{i}<=sTime);
        out = setdiff(1:numel(ClosestValues{i}),in);
        OnIndex = pi{i}(in);
        notCloseOn = pi{i}(out);
        on = [OnIndex',on];
        nc = [notCloseOn',nc];
        
        c=1;
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
        
        c=1;
        for j=1:numel(notCloseOn)
            offset = start-notCloseOn(j);
            offsetE = notCloseOn(j)+after-size(traces,1);
            if offsetE<0
                if offset<0
                    PeakTypeOther(:,c) = traces(notCloseOn(j)-start:notCloseOn(j)+after,i);
                    c=c+1;
                else
                    temp = padarray(traces(:,i),[offset,0],nan,'pre');
                    PeakTypeOther(:,c) = temp(notCloseOn(j)+offset-start+1:notCloseOn(j)+after+offset+1);
                    c=c+1;
                end
                
            else
                temp = padarray(traces(:,i),[offsetE,0],nan,'post');
                PeakTypeOther(:,c) = temp(notCloseOn(j)-start:notCloseOn(j)+after);
                c=c+1;
            end
        end
        
    end
end

end