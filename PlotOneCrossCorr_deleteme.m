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
