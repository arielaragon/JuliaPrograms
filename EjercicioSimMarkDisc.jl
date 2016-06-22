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

P=zeros(4, 4)
P[1, :]=1/4
P[2, :]=[0.0, 0.5, 0.3, 0.2]
P[3, :]=[1/3, 1/3, 1/6, 1/6]
P[4, :]=[0.05, 0.6, 0.25, 0.2]
P
Distributions.sampler
import Distributions.sampler
function sampler(V)
    U=rand()
    X=1
    while cumsum(V)[X]<U
        X=X+1
    end
    X
end
function simulaCM(i, P, n)
    states=zeros(Int, n)
    states[1]=i
    for j in 2:n
        states[j]=sampler(P[states[j-1], :]')
    end
    return states
end
estados=simulaCM(2, P, 1000)
function erase(states, m)
    l=length(estados)
    incomplete=copy(estados)
    e=zeros(Int, m)
    e=round(Int, l*rand(m))
    for i in 1:m
        incomplete[e[i]]=0
    end
    return incomplete
end
incomplete=erase(estados, 100)
function best_t(k, l, d, matriz)
    n=size(matriz)[1]
    V=ones(n^d)
    for i=0:(n^d-1)
        I=zeros(d+2)
        I[1]=k-1
        I[d+2]=l-1
        I[2:(d+1)]=digits(i,n,d)
        I=I+1
        for j=2:(d+2)
            V[i+1]=V[i+1]*matriz[I[j],I[j-1]]
        end
    end
    V=V/sum(V)
    q=sampler(V)  #Activar este renglón hace que me de un muestreo del vector con todas las trayectorias posibles
    #q=indmax(V)   #Activar este renglón hace que me de el estimador más probable
    X_digits=digits(q,n,d)+1
    return X_digits[end-d+1:end]
end
#Funcion que reciba una cadena con datos faltantes (asume que un dato cero es dato faltante) y regresa la cadena completa

function completa(estados, matriz)
    n=size(matriz)[1]
    l=length(estados)
    estados_a=copy(estados)
    aux1=0
    k=1
    if estados_a[1]==0
        while estados_a[k]==0
            aux1=aux1+1
            k=k+1
        end
        maux=matriz^(aux1)
        estados_a[1]=indmax(maux[1:n, estados_a[k]])
        k=1
        aux1=0
    end
    while k<l
        if estados_a[k]!=0
            i=estados_a[k]
            k=k+1
            while estados_a[k]==0
                aux1=aux1+1
                k=k+1
            end
            j=estados_a[k]
            if aux1>0
                estados_a[k-aux1:k-1]=best_t(i, j, aux1, matriz)
            end
            aux1=0
        end
    end
    return estados_a
end
#Solo un ejemplo

b=[1;0;1;0;0;2;1;0;0;0;2]
a=zeros(2, 2)
a[1, 1]=1/3
a[1, 2]=2/3
a[2, 1]=3/4
a[2, 2]=1/4
c=completa(b, a)
p0=mat_disc(incomplete[find(f -> (f>0),incomplete)], 4)
function EME(datos, P0, burn_in, iter)
    d=size(P0)[1]
    n=length(datos)
    m=count(f -> (f>0), datos)
    θ=zeros(d, d, iter)
    θ[:,:,1]=P0
    data_aum=zeros(n)
    for k=2:iter
        data_aum=completa(datos, θ[:,:,k-1])
        θ[:,:,k]=mat_disc(data_aum, d)
    end
    return θ
end
est=EME(incomplete, p0, 1, 1000)

