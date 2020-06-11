covid = csvread('data_branching_process/state_matlab_covid_confirmed.csv',1,0);

[N1 N2]=size(covid);

bins=8;
T=N2;

output1=[];
output2=[];
output3=[];
output4=[];

mus=zeros(N1,1);
alphas=zeros(N1,1);
betas=zeros(N1,1);

K1s=zeros(N1,bins);
Kups=K1s;
Klows=K1s;

for i=1:N1
%for i=2:2
y=covid(i,:);
Nt=y(1:N2);




    
[K1 a1 b1 mu1 p lam AIC Kup Klow]=EM_corona_discrete_hist(Nt,2000,T,bins,200);   

eff_bins=sum(K1>0);
RMSE=mean((lam'-Nt).^2)^.5;
AIC=2*(3+eff_bins)-2*(sum(log(lam'+.000001).*Nt)-sum(lam));
logL=(sum(log(lam'+.000001).*Nt)-sum(lam));



output1=[output1; K1';];
output1=[output1; Kup';];
output1=[output1; Klow';];

output3=[output3;[a1 b1 mu1];];
output4=[output4;[a1 b1 mu1 K1(end) bins logL AIC RMSE];];

output2=[output2; lam';];
output2=[output2; Nt;];

mus(i)=mu1;
alphas(i)=a1;
betas(i)=b1;
K1s(i,:)=K1';
Kups(i,:)=Kup';
Klows(i,:)=Klow';

subplot(1,3,i);
plot([1:max(size(lam))],lam,[1:max(size(lam))],Nt,'o')
drawnow

%     clam=zeros(1,T);
%     cy=zeros(1,T);
%     cy(1)=Nt(1);
%     clam(1)=lam(1);
%     for j=2:T
%     clam(j)=clam(j-1)+lam(j);
%     cy(j)=cy(j-1)+Nt(j);
%     end
%          subplot(2,2,i)
%          plot(clam,cy-clam,'ro',[.01:.01:1]*max(cy),...
%         [.01:.01:1]*max(cy)-[.01:.01:1]*max(cy)...
%         ,'b-',[.01:.01:1]*max(cy),...
%         [.01:.01:1]*max(cy)-ones(100,1)'*1.36*max(cy)^.5-[.01:.01:1]*max(cy),'b-.',...
%         [.01:.01:1]*max(cy),[.01:.01:1]*max(cy)+ones(100,1)'*1.36*max(cy)^.5-[.01:.01:1]*max(cy),'b-.');
%         drawnow;



end

 csvwrite('R0_state_confirmed.csv',output1);
 csvwrite('other_params_state_confirmed.csv',output3);
 csvwrite('intensity_state_confirmed.csv',output2);
 %csvwrite('mod_compare_state_confirmed.csv',output4,'precision',12);
 dlmwrite('mod_compare_state_confirmed.csv', output4, 'delimiter', ',', 'precision', 12);
