function [astro,neuro,a2n]=CalcCorrHist(CosSim, NumNeuro)
[a,~]=(nanmax(CosSim,[],3));
a(logical(a==0))=nan;

astro = triu(a(NumNeuro+1:end,NumNeuro+1:end),1);
astro = astro(:);
astro(astro==0)=[];
astro(isnan(astro))=0;

neuro = triu(a(1:NumNeuro,1:NumNeuro),1);
neuro = neuro(:);
neuro(neuro==0)=[];
neuro(isnan(neuro))=0;

a2n = a(1:NumNeuro,NumNeuro+1:end);
a2n  = a2n(:);
a2n(isnan(a2n))=0;

end