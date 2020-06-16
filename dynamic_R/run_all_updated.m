covid = csvread('matlab_covid_mortality.csv',1,0);

[N1 N2]=size(covid);

bins=9;
T=N2;

output1=[];
output2=[];
output3=[];

mus=zeros(N1,1);
alphas=zeros(N1,1);
betas=zeros(N1,1);
K1s=zeros(N1,bins);
Kups=K1s;
Klows=K1s;


for i=1:N1
y=covid(i,:);
Nt=y(1:N2);


    
[K1 a1 b1 mu1 p lam AIC Kup Klow]=EM_corona_discrete_hist_boundary(Nt,100,T,bins,200);   
   

output1=[output1; K1';];
output1=[output1; Kup';];
output1=[output1; Klow';];

output3=[output3;[a1 b1 mu1];];

output2=[output2; lam';];
output2=[output2; Nt;];

mus(i)=mu1;
alphas(i)=a1;
betas(i)=b1;
K1s(i,:)=K1;
Kups(i,:)=Kup;
Klows(i,:)=Klow;


end

plot(output1')

csvwrite('R0.csv',output1);
csvwrite('other_params.csv',output3);
csvwrite('intensity.csv',output2);

