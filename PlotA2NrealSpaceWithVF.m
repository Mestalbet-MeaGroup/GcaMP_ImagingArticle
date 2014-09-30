% plot a2n links on realspace including n2n of VF elecs 

% [e1wo,e2wo,evalwo,~,~,~,~,~]=ParseConnectionValues(6,100);
[e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValuesWithVF(6,100);
e1=e1(101:end);
e2=e2(101:end);
eval = eval(101:end);
ind = find(eval<min(aval));
e1(ind)=[];
e2(ind)=[];
eval(ind)=[];

PlotResultsOnMEA(e1,e2,eval,mask,ch,a1,a_e2,aval);
