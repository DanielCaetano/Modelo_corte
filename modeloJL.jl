using JuMP
using Cbc

#modelo lendo arquo ntancas

file = ["sinstance_01_1pol_sep.txt","sinstance_01_2pol_sep.txt","sinstance_01_3pol_sep.txt","sinstance_01_4pol_sep.txt","sinstance_01_5pol_sep.txt","sinstance_01_6pol_sep.txt","sinstance_01_7pol_sep.txt","sinstance_01_8pol_sep.txt","sinstance_01_9pol_sep.txt","sinstance_01_10pol_sep.txt"]

function modelo_final(file::Array)
    while file != []
        ele = file[1]
        println(ele)
        entrada = readlines("/home/daniel/Documentos/projetos/julia/resolucao/entradas/$ele")
        vertex = Dict{Int, Any}() # vertices
        edges= Dict{Int, Any}() # arestas
        po=()
        ares=()
        ps =1
        vertex[ps]=po
        edges[1]=ares

        popfirst!(entrada)
        function criar_vertex()
            ps = 1
            println("PS: ",ps)
            for i in entrada
                linha = split(i,' ')
                println(linha)
                p1 = (parse(Float64, linha[1]),parse(Float64, linha[2]))
                p2 = (parse(Float64, linha[3]),parse(Float64, linha[4]))
                t1=true
                t2=true
                for k in 1:length(vertex)
                    #println("p1: $p1 p2: $p2 k:",vertex[k],"\n")
                    if p1 == vertex[k]
                        #println("Repetiu p1:  --------- $p1")
                        t1=false
                    end
                    if p2 == vertex[k]
                        #println("Repetiu p2:  --------- $p2")
                        t2=false
                    end
                end
                if t1==true || t2==true
                    if p1 == p2
                        vertex[ps]=p1
                        ps+=1
                    else
                        if t1 == true
                            vertex[ps]=p1
                            ps+=1
                        end
                        if t2 == true
                            vertex[ps]=p2
                            ps+=1
                        end
                    end
                end
                #println(ponto)
            end

        end

        function criar_edges()
            ps=1
            for i in entrada
                linha = split(i,' ')
                p1 = (parse(Float64, linha[1]),parse(Float64, linha[2]))
                p2 = (parse(Float64, linha[3]),parse(Float64, linha[4]))
                t1=false
                t2=false

                j=1
                ky1=0.0
                ky2=0.0
                while j <= length(vertex)
                    #println("p: $p1 I: $i")
                    if p1 == vertex[j]
                        #println("p1 -------- $j")
                        t1=true
                        ky1=j
                    end
                    if p2 == vertex[j]
                        #println("p2 ------- $j")
                        t2=true
                        ky2=j
                    end

                    if t1 == true && t2 == true
                        #println("$ky1 ------- $ky2")
                        edges[ps]=(ky1, ky2)
                        ps+=1
                        break
                    end
                    j+=1
                end
            end

        end



        criar_vertex()
        #pop!(vertex)
        criar_edges()
        #pop!(edges)

        #Inicio para modelo
        println("MATRIZ ",length(vertex))
        t_cut=1
        t_deslocate =5
        nVertex = length(vertex)
        nedges = length(edges)
        # Matriz para amazenar as distâncias entre os vertex


        dist = zeros(nVertex, nVertex)

        # Definição da posição de cada vértice para o desenho
        posX, posY = [], []
        for i in sort(collect(keys(vertex)))
            posI = vertex[i]
            for j in sort(collect(keys(vertex)))
                posJ = vertex[j]
                dist[i,j]=((posI[1]-posJ[1])^2+(posI[2]-posJ[2])^2 )^(1/2)
            end
            append!(posX, posI[1])
            append!(posY, posI[2])
        end


        #ADJACENTES
        Ac = zeros(nVertex, nVertex)
        # Definição da posição de cada vertice para o desenho

        for i in sort(collect(keys(edges)))
            Ac[edges[i][1], edges[i][2]] = dist[edges[i][1], edges[i][2]]
        end
        U = zeros(nVertex)


        mis = zeros(nVertex, nVertex)

        # Definição da posição de cada vértice para o desenho
        for i in sort(collect(keys(vertex)))
            for j in sort(collect(keys(vertex)))
                mis[i,j]= dist[i,j]/t_deslocate
            end
        end

        #solver Cbc para modelo matematico
        println("MODELO")

        #=
        X -> i,j - todas as edges, t-tempo corte
        Y ->

        dist ij -> distancia entre i e j
        t_cut -> tempo de corte
        t_deslocate -> tempo de movimento

        =#
        # Definição do modelo para a solução, duração maxima de 5 horas
        model=Model(with_optimizer(Cbc.Optimizer, seconds=18000))
        t_inicio = time()
        # Definição da variável matricial Xijy e Yt
        # Note que a variável nVertex foi definida anteriormente após dicionario com as cidades
        @variable(model,X[i=1:nVertex,j=1:nVertex, t=1:2*nedges; i!=j],Bin)
        @variable(model,Y[j=1:nVertex,t=1:2*nedges],Bin)

        # Função objetivo
        @objective(model, Min, sum(X[i, j, t] * dist[i, j]/t_deslocate for i in 1:nVertex, j in 1:nVertex, t in 1:2*nedges if i != j && (Ac[i,j]==0 || Ac[j,i]==0))
                +sum((((sum((X[edges[i][1], edges[i][2], t]+ X[edges[i][2], edges[i][1], t]) for t in 1:2*nedges)-1)
                    *dist[edges[i][1], edges[i][2]])/t_deslocate)+dist[edges[i][1], edges[i][2]]/t_cut for i in 1:length(edges))
                +sum(X[i,j,1] * (((vertex[i][1]-0)^2+(vertex[i][2]-0)^2 )^(1/2)) for i in 1:nVertex for j in 1:nVertex if i != j))

        #(1)
        for t in 1:2*nedges
            @constraint(model, sum(X[i,j,t] for i in 1 : nVertex for j in 1 : nVertex  if i != j ) <= 1)
        end

        #(2)
        for i in 1:nedges
            @constraint(model,sum(X[edges[i][1],edges[i][2],t] + X[edges[i][2],edges[i][1],t]
             for t in 1:2*nedges) >= 1)
        end
        #(3)
        for j in 1:nVertex
            for t in 2:(2*nedges)
                @constraint(model,((sum(X[i,j,t-1] for i in 1:nVertex if i!=j)
                -sum(X[j,k,t] for k in 1:nVertex if k!=j))-Y[j, t]) == 0)
            end
        end

        #Resolução modelo
        JuMP.optimize!(model)
        t_fim = time()
        soma =0
        timeFinal = t_fim - t_inicio
        fo = JuMP.objective_value(model)
        status = JuMP.termination_status(model)

        arq1 = open("/home/daniel/Documentos/projetos/julia/resolucao/resultados/res_$ele", "a")
        write(arq1, "PROBLEM: $ele\n")
        for t in 1:2*nedges
            for i in 1:nVertex
                for j in 1:nVertex
                    if i != j && JuMP.value(X[i, j,t])!=0
                        d=dist[i,j]
                        soma+=d
                        println("Aresta:( $i, $j ), Distance: $d  T: $t ")
                        write(arq1, "A:($i, $j) Distance: $d  T: $t \n")
                    end
                end
            end
        end
        write(arq1, "\nF.O: $fo  ---  Distancia: $soma ---  TIME: $timeFinal  ---  Status: $status\n\n")
        close(arq1)

        #Criar arquivo resultado total
        arq2 = open("/home/daniel/Documentos/projetos/julia/resolucao/resultados/total.txt", "a")
        write(arq2, "PROBLEM: $ele  ---  F.O: $fo  ---  Distancia: $soma ---  TIME: $timeFinal  ---  Status: $status\n")
        close(arq2)

        println("Objective is: $fo")
        println("Time is: ", JuMP.primal_status(model))
        println("Status is: $status")
        println("Tempo total sistema $timeFinal maior que tempo do modelo JuMP: ")
        println("Distancia total: $soma")

        popfirst!(file)
    end
end

modelo_final(file)
