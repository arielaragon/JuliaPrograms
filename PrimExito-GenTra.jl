using Distributions
using StatsBase
e=[0 0 0 1 0 0 ]
#Estima la probabilidad de obtener el primer exito dado un vector, y se le puede indicar como queremos que identifique
#el exito ya sea mediante el 1 o el 0
function PrimExito(e,indext)
  b=length(e)
  c=zeros(length(e))
  v=0
  i=1
  while e[i]==indext
    v=v+1
    i=i+1
  end
  d=((1/b)^(v))*(1-1/b)
end
PrimExito(e,0)
PrimExito(e,1)

vp=[1/5,1/5,3/5]
P=[1/3 1/3 1/3;
   1/2 1/4 1/4;
   0 1/2 1/2]
function camino(numest,vp,Mat,long)
  est=1:numest
  vp=WeightVec(vp)
  path=sample(est,vp)
    for i=2:long
      path=vcat(path,sample(est,WeightVec(vec(Mat[path[end],:]))))
    end
  return path
end
camino(3,vp,P,10)


