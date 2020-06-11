seir<-function(N,T,is,rt,r0,gamma,mu,mort){
  
  #initial cases
  I=numeric(T)
  S=numeric(T)
  R=numeric(T)
  D=numeric(T)
  E=numeric(T)
  dI=numeric(T)
  dD=numeric(T)
  I[1]=is
  S[1]=N-is-rt
  R[1]=rt
  E[1]=0
  NewCase=numeric(T)
  
  beta=(r0*gamma) # 
  
  
  # IR simulation
  for(i in 2:T){
    
    I[i]=I[i-1]+mu*E[i-1]-gamma*I[i-1]
    dI[i]=mu*E[i-1]-gamma*I[i-1]
    dD[i]=mort*gamma*I[i-1]
    R[i]=R[i-1]+gamma*I[i-1]
    E[i]=E[i-1]+beta*I[i-1]*S[i-1]/N-mu*E[i-1]
    S[i]=S[i-1]-beta*I[i-1]*S[i-1]/N
    D[i]=mort*R[i]
    # new case each day
    NewCase[i]=mu*E[i-1]
  }
  output=c()
  output$dI=dI
  output$dD=dD
  output$I=I
  output$muE=mu*E
  return(output)
  
}
