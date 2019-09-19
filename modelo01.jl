using JuMP

#modelo lendo arquo ntancas

entrada = readlines("/home/daniel/Documentos/projetos/julia/resolucao/entradas/sinstance_01_3pol_sep.txt")


pontos = Dict{Int, Any}()
arestas= Dict{Int, Any}()
passos = Dict{Any, Any}()
po=()
ares=()
ps =1
pontos[ps]=po
arestas[1]=ares
passos[1]=()

#= #definir string para poder utilizar split(String , Char)
whole::String
whole="this^DB+Verb+Pos^DB+Noun+Inf2+A3sg+Pnon+Nom"
println(split(whole,['^','+']))
=#
#campos=split(entrada,' ')
#println(campos)

popfirst!(entrada)
function criar_pontos()
    global arestas
    global ps = 1
    for i in entrada
        linha = split(i,' ')
        println(linha)
        p1 = (parse(Float64, linha[1]),parse(Float64, linha[2]))
        p2 = (parse(Float64, linha[3]),parse(Float64, linha[4]))
        t1=true
        t2=true
        for k in 1:length(pontos)
            #println("p1: $p1 p2: $p2 k:",pontos[k],"\n")
            if p1 == pontos[k]
                #println("Repetiu p1:  --------- $p1")
                t1=false
            end
            if p2 == pontos[k]
                #println("Repetiu p2:  --------- $p2")
                t2=false
            end
        end
        if t1==true || t2==true
            if p1 == p2
                pontos[ps]=p1
                ps+=1
            else
                if t1 == true
                    pontos[ps]=p1
                    ps+=1
                end
                if t2 == true
                    pontos[ps]=p2
                    ps+=1
                end
            end
        end
        #println(ponto)
    end

end

function criar_arestas()
    global arestas
    global ps=1
    for i in entrada
        linha = split(i,' ')
        p1 = (parse(Float64, linha[1]),parse(Float64, linha[2]))
        p2 = (parse(Float64, linha[3]),parse(Float64, linha[4]))
        t1=false
        t2=false

        j=1
        ky1=0.0
        ky2=0.0
        while j <= length(pontos)
            #println("p: $p1 I: $i")
            if p1 == pontos[j]
                #println("p1 -------- $j")
                t1=true
                ky1=j
            end
            if p2 == pontos[j]
                #println("p2 ------- $j")
                t2=true
                ky2=j
            end

            if t1 == true && t2 == true
                #println("$ky1 ------- $ky2")
                arestas[ps]=(ky1, ky2)
                ps+=1
                break
            end
            j+=1
        end
    end

end



criar_pontos()
#pop!(pontos)
criar_arestas()
#pop!(arestas)

#= #Criar pontos no arquivo
arq = open("/home/daniel/Documentos/projetos/julia/resolucao/resposta.txt", "a")
println("Pontos")
x=1

#write(arq, "$x")
#close(arq)

for i in 1:length(pontos)
    t1=pontos[i][1]
    t2=pontos[i][2]
    write(arq, "$t1, $t2\n")
    println(pontos[i])
end

println("Arestas")
for i in 1:length(arestas)
    t1=arestas[i][1]
    t2=arestas[i][2]
    write(arq, "$t1, $t2\n")
    #println(i)
end

close(arq)
=#

#Inicio para modelo
println("MATRIZ")
t_corte=1
t_movimento =5
nPontos = length(pontos)
nArestas = length(arestas)
# Matriz para amazenar as distâncias entre os pontos


dist = zeros(nPontos, nPontos)

# Definição da posição de cada vértice para o desenho
posX, posY = [], []
for i in sort(collect(keys(pontos)))
    posI = pontos[i]
    for j in sort(collect(keys(pontos)))
        posJ = pontos[j]
        dist[i,j]=((posI[1]-posJ[1])^2+(posI[2]-posJ[2])^2 )^(1/2)
    end
    append!(posX, posI[1])
    append!(posY, posI[2])
end


#ADJACENTES
Ac = zeros(nPontos, nPontos)
# Definição da posição de cada vertice para o desenho

for i in sort(collect(keys(arestas)))
    Ac[arestas[i][1], arestas[i][2]] = dist[arestas[i][1], arestas[i][2]]
end
U = zeros(nPontos)


move = zeros(nPontos, nPontos)

# Definição da posição de cada vértice para o desenho
for i in sort(collect(keys(pontos)))
    for j in sort(collect(keys(pontos)))
        move[i,j]= dist[i,j]/t_movimento
    end
end
#println(U)
#solver Cbc para modelo matematico
using Cbc
println("MODELO")



#=
dist ij -> distancia entre i e j
t_corte -> tempo de corte
t_movimento ->

=#
# Definição do modelo para a solução
model=Model(with_optimizer(Cbc.Optimizer, seconds=18000))

# Definição da variável matricial xij
# Note que a variável nPontos foi definida anteriormente após dicionario com as cidades
@variable(model,x[i=1:nPontos,j=1:nPontos, t=1:2*nArestas; i!=j],Bin)
@variable(model,Y[j=1:nPontos,t=1:2*nArestas],Bin)

# Função objetivo
@objective(model, Min, sum(x[i, j, t] * dist[i, j]/t_movimento for i in 1:nPontos, j in 1:nPontos, t in 1:2*nArestas if i != j)
        +(sum(Ac[arestas[i][1], arestas[i][2]]/t_corte - Ac[arestas[i][1], arestas[i][2]]/t_movimento for i in 1:length(arestas)))/2
        +sum(x[i, j, 1] * (((pontos[i][1]-0)^2+(pontos[i][2]-0)^2 )^(1/2)) for i in 1:nPontos, j in 1:nPontos if i != j))

#(0)
#@constraint(model, sum(x[1,j,t] for i in 1 : nPontos for j in 1 : nPontos  if i != j ) <= 1)

#(1)
for t in 1:2*nArestas
    @constraint(model, sum(x[i,j,t] for i in 1 : nPontos for j in 1 : nPontos  if i != j ) <= 1)
end

#(2)
for i in 1:nArestas
    @constraint(model,sum(x[arestas[i][1],arestas[i][2],t] + x[arestas[i][2],arestas[i][1],t]
     for t in 1:2*nArestas) >= 1)
end
#(3)
for j in 1:nPontos
    for t in 2:(2*nArestas)
        @constraint(model,((sum(x[i,j,t-1] for i in 1:nPontos if i!=j)
        -sum(x[j,k,t] for k in 1:nPontos if k!=j))-Y[j, t]) == 0)
    end
end


#status = solve(model)

#Resolução modelo
JuMP.optimize!(model)

#if status == :Optimal
edgeOrigin = []
edgeDest = []
#println(su)
soma =0
for t in 1:2*nArestas
    global soma
    for i in 1:nPontos
        for j in 1:nPontos
            if i != j && JuMP.value(x[i, j,t])!=0

                append!(edgeOrigin, i)
                append!(edgeDest, j)
                d=dist[i,j]
                soma+=d
                println("Aresta:( $i, $j ), Distance: $d  T: $t ")
            end
        end
    end
end
println("Distancia total: $soma")
