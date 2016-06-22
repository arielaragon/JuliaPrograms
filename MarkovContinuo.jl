using StatsBase
using Distributions
#example of trajectory
tray0=[0.1 2;
      0.15 3;
      0.25 1;
      0.40 2;
      0.67 1;
      0.88 4;
      1.2 2]
# infinitesimal generator matrix
Q=[-4 2 1 1;
   2 -6 3 1;
   1 1 -4 2;
   1 3 2 -6]

#Los siguentes dos programitas estiman, primero; el numnero de saltos de un estado a otro, un vector para indicar cuantas veces
#aparece el estado, con los puentes se le puede aplicar
function cuenta(trayecto)
  tra=trayecto
  tam=size(tra)[1]
  E=tra[1:end,2]
  n=size(unique(E),1)
  T=tra[1:end,1]
  R=zeros(n)
  N=zeros(n,n)
  for i in 1:(tam-1)
    delta=T[i+1]-T[i]
    R[E[i]]=R[E[i]]+delta
    N[E[i],E[i+1]]=N[E[i],E[i+1]]+1
  end
  return (R,N)
end

function EstimMCC(trayecto)
  tra= trayecto
  tam=size(tra)[1]
  E=tra[1:end,2]
  n=size(unique(E),1)
  T=tra[1:end,1]
  R=zeros(n)
  N=zeros(n,n)
  Q=zeros(n,n)
  P=zeros(n)
  for i in 1:(tam-1)
    delta=T[i+1]-T[i]
    R[E[i]]=R[E[i]]+delta
    N[E[i],E[i+1]]=N[E[i],E[i+1]]+1
  end
  Q=N./R
    for i in 1:n
        Q[i,i]=-sum(Q,2)[i]
    end
    P=N./sum(N,2)
  return (Q,P)
end





function MJPBisect(ini,fin,m,E,Q,Qii,T1,T2,cond)
  tray=[T1 ini;T2 fin]
  t0=T2-T1
  t=t0/2
  t1=-t
  # Transition matrix of time lag t
  P=expm(t*Q)
  # We assign the probabilities of transitioning through each state
  prob=zeros(m)
  for i=E
    prob[i]=P[ini,i]*P[i,fin]
  end

  # We modify the probabilities if we are being conditioned to have at least 2 jumps
  if cond
    if ini==fin
      prob[ini]-=exp(-t0*Qii[ini])
    elseif Qii[ini]==Qii[fin]
      r=exp(t1*Qii[ini])*t*exp(t1*Qii[ini])*Q[ini,fin]
      prob[ini]-=r
      prob[fin]-=r
    else
      r=exp(t1*Qii[ini])
      r=r*(r-exp(t1*Qii[fin]))*Q[ini,fin]/(Qii[fin]-Qii[ini])
      prob[ini]-=r
      prob[fin]-=r
    end
    # We sample to see wich state we passed through
    mid=sample(E,weights(prob))

    # Once we have the state, we check to see in which manner we did so (number of jumps)
    if ini==mid
      if mid==fin
        r=exp(t1*Qii[ini])
        rc=P[ini,mid]-r
        (bef,aft)=sample([(0,2),(2,0),(2,2)],weights([r*rc,rc*r,rc*rc]))
      else
        r=exp(-t*Qii[ini])
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(0,2),(2,1),(2,2)],weights([r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    else
      if mid==fin
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        r1=exp(t1*Qii[mid])
        (bef,aft)=sample([(1,2),(2,0),(2,2)],weights([r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      else
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(1,1),(1,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    end
  else
    # We sample to see wich state we passed through
    mid=sample(E,weights(prob))

    # Once we have the state, we check to see in which manner we did so (number of jumps)
    if ini==mid
      if mid==fin
        r=exp(t1*Qii[ini])
        rc=P[ini,mid]-r
        (bef,aft)=sample([(0,0),(0,2),(2,0),(2,2)],weights([r*r,r*rc,rc*r,rc*rc]))
      else
        r=exp(t1*Qii[ini])
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(0,1),(0,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    else
      if mid==fin
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        r1=exp(t1*Qii[mid])
        (bef,aft)=sample([(1,0),(1,2),(2,0),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      else
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(1,1),(1,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    end
  end

  if bef==1
    if Qii[fin]>Qii[ini]
      tray=vcat(tray,[T1+rand(Truncated(Exponential(1/(Qii[fin]-Qii[ini])),0,t)) mid])
    elseif Qii[fin]<Qii[ini]
      tray=vcat(tray,[T1+t-rand(Truncated(Exponential(1/(Qii[ini]-Qii[fin])),0,t)) mid])
    else
      tray=vcat(tray,[T1+rand(1)*t mid])
    end
  elseif bef==2
    tray=vcat(MJPBisect(ini,mid,m,E,Q,Qii,T1,T1+t,true),tray)
  end
  if aft==1
    if Qii[fin]>Qii[ini]
      tray=vcat(tray,[T1+t+rand(Truncated(Exponential(1/(Qii[fin]-Qii[ini])),0,t)) fin])
    elseif Qii[fin]<Qii[ini]
      tray=vcat(tray,[T2-rand(Truncated(Exponential(1/(Qii[ini]-Qii[fin])),0,t)) fin])
    else
      tray=vcat(tray,[T1+t+rand(1)*t fin])
    end
  elseif aft==2
    tray=vcat(tray,MJPBisect(mid,fin,m,E,Q,Qii,T1+t,T2,true))
  end
  return(tray)
end
#####################################################
function MJPBridgeBisection(tray,Q)
  m=length(Q[1,1:end])
  n=length(tray[1:end,1])-1
  E=1:m
  Qii=-diag(Q)
  tray0=hcat([],[])
  for i = 1:n
    mat=MJPBisect(tray[i,2],tray[i+1,2],m,E,Q,Qii,tray[i,1],tray[i+1,1],false)
    mat=sortrows(mat)
    rng=hcat(mat[1,1],mat[1,2])
    for k = 2:length(mat[1:end,1])
      if mat[k,2]!=mat[k-1,2]
        rng=vcat(rng,mat[k,1:2])
      end
    end
    tray0=vcat(tray0,rng)
  end
  tray0=vcat(tray0,tray[end;1:2])
  return(tray0)
end

@time MJPBridgeBisection(tray0,Q)


