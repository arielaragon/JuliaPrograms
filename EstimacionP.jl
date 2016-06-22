#Dado un vector de estados, el programita identifica el numero de estados, y estima la matriz P
function mat_disc(estados)
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
        m[i, 1:n]=m[i, 1:n]/s[i]
    end
    return m
end

estados=[1, 1, 3, 2, 4, 2, 3, 1, 2, 4 ,2, 3 , 1, 3, 2, 4, 3, 2, 3, 3, 2, 1, 2, 1, 4, 2, 3, 4, 2, 2, 3, 4,3, 4, 2, 1, 2, 2, 1, 1 ,3, 2, 3, 4, 3]
m=mat_disc(estados)

