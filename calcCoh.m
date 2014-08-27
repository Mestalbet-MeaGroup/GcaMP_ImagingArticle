function [freq,coh,UpperSig,LowerSig]=calcCoh(source1,source2,fs)

leng = min([numel(source1),numel(source2)]);
freq=((1:leng)-1)/(leng-1)*fs/2;
% fftsource1=pwelch(source1);
% fftsource2=pwelch(source2);
% fftsource1=fftsource1(1:leng);
% fftsource2=fftsource2(1:leng);
        
%% Calc Spectrum and Cross-Spectrum

[cs,freq] = cpsd(source1,source2,[],[],[],fs);
[spectSource1,~] = cpsd(source1,source1,[],[],[],fs);
[spectSource2,~] = cpsd(source2,source2,[],[],[],fs);
% spectSource1=fftsource1.*conj(fftsource1);
% spectSource2=fftsource2.*conj(fftsource2);

%% Calculate Coherence and signficance values
coh=abs(cs)./sqrt(spectSource1.*spectSource2); % No smoothing
nu=numel(freq);
UpperSig=2*sqrt(ones(leng,1)*1/nu.*(1./coh(:).^2 - 1));
LowerSig=-2*sqrt(ones(leng,1).*1/nu.*(1./coh(:).^2 - 1));

% nu= 2./sum(smooth.^2);  % 2*number of epochs*1./sum(smooth.^2);
% s=ones(leng,1)*sqrt(1-(0.05)^(2/(nu-2))); % 0.05% significance
% psu=2*sqrt(ones(leng,1)*1/nu.*(1./coh(:).^2 - 1));
% psl=-2*sqrt(ones(leng,1).*1/nu.*(1./coh(:).^2 - 1));
% freq=freq';

end
