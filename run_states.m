covid = csvread('data_branching_process/state_matlab_covid_mortality.csv',1,0);

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

y=covid(i,:);
Nt=y(1:N2);


    
 [K1 a1 b1 mu1 p lam AIC Kup Klow]=EM_corona_discrete_hist(Nt,1000,T,bins,100);   


eff_bins=sum(K1>0);
AIC=2*(3+eff_bins)-2*(sum(log(lam'+.000001).*Nt)-sum(lam));
RMSE=mean((lam'-Nt).^2)^.5;
logL=(sum(log(lam'+.000001).*Nt)-sum(lam));

output1=[output1; K1';];
output1=[output1; Kup';];
output1=[output1; Klow';];

output3=[output3;[a1 b1 mu1];];
output4=[output4;[a1 b1 mu1 K1(end) logL AIC RMSE];];

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


end

csvwrite('R0_state.csv',output1);
csvwrite('other_params_state.csv',output3);
csvwrite('intensity_state.csv',output2);
csvwrite('mod_compare_state.csv',output4);

