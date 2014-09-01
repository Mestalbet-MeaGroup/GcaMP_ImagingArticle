% for i=1:9
%     loadcult(i);
%     [atrigpsth{i},atriglag{i}]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);
%     [ntrigpsth{i},ntriglag{i}]=CalcPSTHneurotriggers(t,ic,traces,time,bs,be);
%     clear_all_but('atrigpsth','atriglag','ntrigpsth','ntriglag','i');
% end

load('PSTHdata.mat')
subplot(2,1,1)
hold all;
for i=1:9
   plot(atriglag{i},nansum(zscore(atrigpsth{i},0,2),1)); 
end
title('Astrocyte Triggered PSTH');
subplot(2,1,2)
hold all;
for i=1:9
   plot(ntriglag{i},nansum(zscore(ntrigpsth{i},0,2),1)); 
end
title('Burst Triggered Astrocyte Peak PSTH');