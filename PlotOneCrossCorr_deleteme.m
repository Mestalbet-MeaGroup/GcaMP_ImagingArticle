subplot(4,1,1)
plotyy(DataSet{6}.dfTime,zscore(DataSet{6}.GFR),DataSet{6}.dfTime,zscore(DataSet{6}.dfTraces(:,80)));
subplot(4,1,2)
plotyy(DataSet{6}.dfTime,zscore(filtfilt(MakeGaussian(0,30,120),1,DataSet{6}.GFR)),DataSet{6}.dfTime,zscore(filtfilt(MakeGaussian(0,30,120),1,DataSet{6}.dfTraces(:,80))) )
subplot(4,1,3)
[a,b]=xcorr(zscore(DataSet{6}.GFR),zscore(DataSet{6}.dfTraces(:,80)),'unbiased');
tgfr = xcorr(zscore(DataSet{6}.GFR),zscore(DataSet{6}.GFR),'unbiased');
ttr = xcorr(zscore(DataSet{6}.dfTraces(:,80)),zscore(DataSet{6}.dfTraces(:,80)),'unbiased');
plot(b,a);
axis tight;
xlim([-1,1]*10^4);
subplot(4,1,4)
[a,b]=xcorr(zscore(filtfilt(MakeGaussian(0,30,120),1,DataSet{6}.GFR)),zscore(filtfilt(MakeGaussian(0,30,120),1,DataSet{6}.dfTraces(:,80))),'unbiased');
plot(b,a);
axis tight;
xlim([-1,1]*10^4);


%% Plot Corr with different normalizations
w=1;
wt=300-238;
% w=6;
% wt=80;
subplot(5,1,1)
[a1,b1]=xcorr(zscore(DataSet{w}.GFR),zscore(DataSet{w}.dfTraces(:,wt)),'unbiased');
plot(b1,a1);
axis tight;
title('Unbiased Zscore');
subplot(5,1,2)
[a2,b2]=xcorr(DataSet{w}.GFR,DataSet{w}.dfTraces(:,wt),'unbiased');
plot(b2,a2);
axis tight;
title('Unbiased');
subplot(5,1,3)
[a3,b3]=xcorr(DataSet{w}.GFR,DataSet{w}.dfTraces(:,wt),'coef');
plot(b3,a3);
axis tight;
title('Coef');
subplot(5,1,4)
[a4,b4]=xcorr(zscore(DataSet{w}.GFR),zscore(DataSet{w}.dfTraces(:,wt)),'coef');
plot(b4,a4);
axis tight;
title('zscore coef');
subplot(5,1,5)
[b5,a5] = CalcCorrRawFull(zscore(DataSet{w}.GFR),zscore(DataSet{w}.dfTraces(:,wt)),numel(DataSet{6}.GFR));
plot(b5,a5);
axis tight;
title('Unbiased Normalized');

