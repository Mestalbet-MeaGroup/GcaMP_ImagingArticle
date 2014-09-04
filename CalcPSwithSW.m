function ps = CalcPSwithSW(t,ic,step,window);
%Calculates the phase stability with a sliding window

winstart = min(t);
numc = size(ic,2);
c=1;
ps=zeros(numc,numc,numel(min(t):step*12000:max(t)));
while winstart<(max(t)-step*12000)
    winstart = winstart+step*12000;
    [tNew,icNew]=CutSortChannel2(t,ic,winstart,winstart+window*12000);
    temp=CalcPhaseStability(tNew,icNew); 
    locs = ismember(ic(1,:),icNew(1,:));
    ps(locs,locs,c)=temp;
%     if size(temp,1)<numc
%         ps(:,:,c) = padarray(temp,[numc-size(temp,1),numc-size(temp,2)],0,'post');
%     else
%         ps(:,:,c) = temp;
%     end
    c=c+1;
end

end