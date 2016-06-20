using StatsBase
using Distributions
function mat_disc(estados, n)
    m=zeros(Float64, n, n)
    l=length(estados)
    for i in 1:l-1
        m[estados[i], estados[i+1]]=m[estados[i], estados[i+1]]+1
    end
    s=zeros(Float64, n)
    for i in 1:n
        for j in 1:n
            s[i]=s[i]+m[i, j]
        end
        m[i, 1:n]=m[i, 1:n]/s[i]
    end
    return m
end
estados=[1, 1, 3, 2, 4, 2, 3, 1, 2, 4 ,2, 3 , 1, 3, 2, 4, 3, 2, 3, 3, 2, 1, 2, 1, 4, 2, 3, 4, 2, 2, 3, 4,3, 4, 2, 1, 2, 2, 1, 1 ,3, 2, 3, 4, 3]
m=mat_disc(estados, 4)
function EME(datos, P0, burn_in, iter)
    d=size(P0)[1]
    n=length(datos)
    m=count(f -> (f>0), datos)
    θ=zeros(d, d, iter)
    θ[:,:,1]=P0
    data_aum=zeros(n)
    for k=2:iter
        θ[:,:,k]=mat_disc(data_aum, d)
    end
    return θ
end

P=zeros(4, 4)
P[1, :]=1/4
P[2, :]=[0.0, 0.5, 0.3, 0.2]
P[3, :]=[1/3, 1/3, 1/6, 1/6]
P[4, :]=[0.05, 0.6, 0.25, 0.2]
P
da=data([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; 2 3 3 1 2 1 2 1 2 2 1 1 3 4 4 1 4 3 2 1])'
est=da[:,2]
##########################
#########MAL##############
EME(est,P,1,10)###########
##########################
##########################

###Con esto se puede hacer una funcion que te genere matrices de estados, que lo que deberia de ser
###es una matriz a la que le pases n vectores aleatorios de estados y luego puedas sacar la matriz
###de transicion.
iter=6
d=size(P)[1]
n=length(est)
θ=zeros(d, d, iter)
θ[:,:,1]=P
data_aum=zeros(n)
for k=2:iter
    θ[:,:,k]=θ[:,:,k-1]+mat_disc(est,d)  ####rand(1:d, d, d)
end
return θ
t= P
while k<=iter
  t[:,:,k]=mean(t[:,:,k-1])
  k=k+1
end

c=3
i=3
B=ones(3,3)
A=zeros(c,c,i)
A[:,:,1]=B
for k=2:i
  A[:,:,k]=A[:,:,k-1]+1
end
return A

#################################################################################################
#################################################################################################
#################################################################################################
da=data([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; 2 3 3 1 2 1 2 5 2 2 1 1 3 4 4 5 4 3 2 5])'


est=da[:,2]
size(unique(est),1)
function Clasif(estados)
    n=size(unique(estados),1)
    m=zeros(Float64, n, n)
    l=length(estados)
    for i in 1:l-1
        m[estados[i], estados[i+1]]=m[estados[i], estados[i+1]]+1
    end
    s=zeros(Float64, n)
    for i in 1:n
        for j in 1:n
            s[i]=s[i]+m[i, j]
        end
        m[i, 1:n]=m[i, 1:n]
    end
    return m
end
k=Clasif(est)
n=5
y=zeros(Float64,n,n)
for i in 1:n
  for j in 1:n
    y[i,j]=mean(k[i,j])
  end
end
return y


n=5
t=zeros(Float64, n)
for i in 1:n
    for j in 1:n

    k[i, 1:n]=k[i, 1:n]
        end
       t[i]=mean(k[i, 1:n])
    end
return t

F=zeros(Float64,n,n)
for i in i:n
  for j in 1:n
     F[i, 1:n]=F[i,1:n]+m[i, 1:n]/s[i]
  end
  return F
end

function mat_disc(estados, n)
    m=zeros(Float64, n, n)
    l=length(estados)
    for i in 1:l-1
        m[estados[i], estados[i+1]]=m[estados[i], estados[i+1]]+1
    end
    s=zeros(Float64, n)
    for i in 1:n
        for j in 1:n
            s[i]=s[i]+m[i, j]
        end
        m[i, 1:n]=(m[i, 1:n]/s[i])/t[i]
    end
    return m
end

mat_disc(est,5)
G=zeros(Float64,n,n)
for i in 1:n
  for j in 1:n
    t[i]=t[i]+G[i,j]
  end
  G[i,j]=G[i,j]+m[i,j]/t[i]
  return G
end


