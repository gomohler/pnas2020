function [K0 alpha beta mu p lam AIC Kup Klow]=EM_corona_discrete_hist_boundary(t,emiter,T,bins,Nconf)

N=max(size(t));

% the end time boundary biases estimates (because we can't observe future
% offspring.  so only calculate p_ij when j is extra days before the end
% time T or earlier
extra=7;

% p is a matrix storing branching probabilities
p=zeros(N,N);

%initial guesses for parameters
K0=.5*ones(bins,1);

alpha=6;
beta=3;
mu=10;


% number of EM iterations
for k=1:emiter


% E-Step
[p lam]=updatep(t,p,K0,alpha,beta,mu,T,bins,extra);   


% M-Step
[K0 alpha beta mu]=updatepar(K0,alpha,beta,t,p,T,bins,extra);

%stabilize weibull
if(beta>4)
    beta=4;
end
% if(mu>1)
%     mu=1;
% end

if(mod(k,10)==0)
   [k/emiter mu alpha beta K0(end)]
end

end

[Kup Klow]=K_confidence(K0,alpha,beta,t,p,T,bins,Nconf);

AIC=2*(bins+3)-2*(sum(log(lam'+.000001).*t)-sum(lam));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate branching probabilities

function [p lam]=updatep(t,p,K0,alpha,beta,mu,T,bins,extra)

N=max(size(t));

lam=zeros(N,1);

for i=1:N
    for j=1:(i-1)
        
            d=i-j;            
            % probability i triggered by j is proportional to triggering
            % kernel evaluated at inter-point times and distances
            p(i,j)=K0(ceil(bins*j/T))*wblpdf(d,alpha,beta)*t(j);
         %   p(i,j)=K0(ceil(bins*j/T))*wbl_new(i,j,alpha,beta)*t(j);


    end
    
    %probablity i is background event proportional to mu background rate
    p(i,i)=mu;
    
    % save intensity at each event for analysis 
    lam(i)=sum(p(i,1:i));
    
    %normalize probabilities to sum to 1
    p(i,1:i)=p(i,1:i)/sum(p(i,1:i));
    

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate parameters from branching probabilities

function [K0 alpha beta mu]=updatepar(K0,alpha,beta,t,p,T,bins,extra)

N=max(size(t));

mu=0;

K0=zeros(size(K0));
Nc=zeros(size(K0));

time_sample=[];
p_sample=[];

for i=1:N

    for j=1:min((i-1),N-extra)

        % parameters are determined by weighted sample mean
        % of inter-point times and square distances        
        
        time_sample=[time_sample; (i-j);]; 
        p_sample=[p_sample; p(i,j)*t(i);];
      
    end
    
    for j=1:(i-1)
        
        K0(ceil(bins*j/T))=K0(ceil(bins*j/T))+p(i,j)*t(i);
        
        
    end
    
    mu=mu+p(i,i)*t(i);
  
end

for i=1:N
  %this is the poor man's version, not integrating over the day that is observation i  
  Nc(ceil(bins*i/T))=Nc(ceil(bins*i/T))+t(i)*(1-exp(-((T-i)/alpha)^beta));
end





K0=K0./(Nc+.000001);

[coef,~] = wblfit(time_sample,[],[],p_sample);
alpha=coef(1);
beta=coef(2);

% disp(string(alpha))
% disp(string(beta))
% K0
% pause(2)


% we don't have a good estimate of K0 at the boundary so just set it
% equal to K0 at previous bin
%K0(end)=K0(end-1);

mu=mu/T;

end

function [Kup Klow]=K_confidence(K0,alpha,beta,t,p,T,bins,Nconf)

Kconf=zeros(bins,Nconf);

for nc=1:Nconf
  
    K0=zeros(size(K0));
    Nc=zeros(size(K0));

    for i=1:T
    for j=1:(i-1)
        
        K0(ceil(bins*j/T))=K0(ceil(bins*j/T))+binornd(t(i),p(i,j));
        
        
    end
   
    Nc(ceil(bins*i/T))=Nc(ceil(bins*i/T))+t(i)*(1-exp(-((T-i)/alpha)^beta));
  
    end

    K0=K0./(Nc+.000001);
    
    % we don't have a good estimate of K0 at the boundary so just set it
    % equal to K0 at previous bin
   %K0(end)=K0(end-1);
    Kconf(:,nc)=K0;
    

end

Kup=zeros(bins,1);
Klow=zeros(bins,1);

for j=1:bins
    Kup(j)=quantile(Kconf(j,:),.975);
    Klow(j)=quantile(Kconf(j,:),.025);
    
end

end


function p=pois(S)

if(S<=100)
temp=-S;
L=exp(temp);
k=0;
p=1; 
while(p > L)
k=k+1;
p=p*rand();
end
p=k-1;
else
p=floor(S+S^.5*randn());
end
end