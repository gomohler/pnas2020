sir<-function(N,T,is,rt,r0,gamma,mort){
  
  #initial cases
  I=numeric(T)
  S=numeric(T)
  R=numeric(T)
  D=numeric(T)
  dI=numeric(T)
  dD=numeric(T)
  I[1]=is 
  S[1]=N-is-rt
  R[1]=rt
  NewCase=numeric(T)
  
  beta=(r0*gamma) # 
  
  
  # IR simulation
  for(i in 2:T){
    
    I[i]=I[i-1]+beta*S[i-1]*I[i-1]/N-gamma*I[i-1]
    dI[i]=beta*S[i-1]*I[i-1]/N
    dD[i]=mort*gamma*I[i-1]
    R[i]=R[i-1]+gamma*I[i-1]
    S[i]=S[i-1]-beta*I[i-1]*S[i-1]/N
    D[i]=mort*R[i]
    # new case each day
    NewCase[i]=beta*I[i-1]*S[i-1]/N
  }
  output=c()
  output$dI=dI
  output$dD=dD
  output$I=I
  return(output)
  
}
