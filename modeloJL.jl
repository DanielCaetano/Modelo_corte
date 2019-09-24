using JuMP
using Cbc

#modelo lendo arquo ntancas
lista = ["rinstance_01_2pol.txt","rinstance_01_3pol.txt","sinstance_01_2pol_sep.txt","sinstance_01_3pol_sep.txt"]

lista1 = ["rinstance_01_4pol.txt","rinstance_01_5pol.txt","sinstance_01_5pol_sep.txt"]

lista2 = ["sinstance_01_2pol_sep.txt","sinstance_01_3pol_sep.txt","sinstance_01_4pol_sep.txt"]

function modelo_final(file::Array)
    while file != []
        ele = file[1]
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
            #println("PS: ",ps)
            for i in entrada
                linha = split(i,' ')
                println(linha)
                p1 = (parse(Float64, linha[1]),parse(Float64, linha[2]))
                p2 = (parse(Float64, linha[3]),parse(Float64, linha[4]))
                t1=true
                t2=true
                for k in 1:length(vertex)
                    if p1 == vertex[k]
                        t1=false
                    end
                    if p2 == vertex[k]
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
                    if p1 == vertex[j]
                        t1=true
                        ky1=j
                    end
                    if p2 == vertex[j]
                        t2=true
                        ky2=j
                    end

                    if t1 == true && t2 == true
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
        println("MATRIZ ")
        t_cut=1
        t_deslocate =5
        nVertex = length(vertex)
        nEdges = length(edges)
        # Matriz para amazenar as distâncias entre os vertex
        dist = zeros(nVertex, nVertex)
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

        #Arestas de corte
        #Eges of cut
        Ac = zeros(nVertex, nVertex)
        # Definição da posição de cada vertice para o desenho
        # Setting the position of each vertex for the drawing

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

        pis = zeros(nEdges, nEdges)
        # Definição da posição de cada vértice para o desenho
        for i in sort(collect(keys(edges)))
                pis[edges[i][1], edges[i][2]] = dist[edges[i][1], edges[i][2]]/t_cut
                #println(pis[edges[i][1], edges[i][2]])
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
        @variable(model,X[i=1:nVertex,j=1:nVertex, t=1:2*nEdges; i!=j],Bin)
        @variable(model,Y[j=1:nVertex,t=1:2*nEdges],Bin)

        # Função objetivo
        @objective(model, Min, sum(X[i, j, t] * mis[i, j] for i in 1:nVertex for j in 1:nVertex for t in 1:2*nEdges if i != j && (Ac[i,j]==0 && Ac[j,i]==0))
                +sum(((sum((X[edges[i][1], edges[i][2], t]+ X[edges[i][2], edges[i][1], t]) for t in 1:2*nEdges)-1)
                    *mis[edges[i][1], edges[i][2]])+pis[edges[i][1], edges[i][2]] for i in 1:nEdges)
                +sum(X[i,j,1] * (((vertex[i][1]-0)^2+(vertex[i][2]-0)^2 )^(1/2)) for i in 1:nVertex for j in 1:nVertex if i != j))

        #(1)
        for t in 1:2*nEdges
            @constraint(model, sum(X[i,j,t] for i in 1 : nVertex for j in 1 : nVertex  if i != j ) <= 1)
        end

        #(2)
        for i in 1:nEdges
            @constraint(model,sum(X[edges[i][1],edges[i][2],t] + X[edges[i][2],edges[i][1],t]
             for t in 1:2*nEdges) >= 1)
        end
        #(3)
        for j in 1:nVertex
            for t in 2:(2*nEdges)
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
        for t in 1:2*nEdges
            for i in 1:nVertex
                for j in 1:nVertex
                    if i != j && JuMP.value(X[i, j,t])!=0
                        d=dist[i,j]
                        soma+=d
                        v_mis = mis[i,j]
                        println("Aresta:( $i, $j ), Distance: $d  T: $t mis: $v_mis")
                        write(arq1, "A:($i, $j) Distance: $d  T: $t mis: $v_mis\n")
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

        fo1 = 0
        fo2 = 0
        fo3 = 0
        for i in 1:nVertex
            for j in 1:nVertex
                for t in 1:2*nEdges
                    if i != j && (Ac[i,j]==0 && Ac[j,i]==0)
                        if JuMP.value(X[i, j, t] ) != 0.0
                            println("X($i, $j, $t)")
                        end
                        fo1 = fo1 + (JuMP.value(X[i, j, t])* mis[i, j])
                    end
                end
            end
        end

        fo2_s = 0
        for i in 1:nEdges
            fo2_s=0
            for t in 1:2*nEdges
                fo2_s =fo2_s+ JuMP.value(X[edges[i][1], edges[i][2], t])+ JuMP.value(X[edges[i][2], edges[i][1], t])
            end
            a1=edges[i][1]
            a2=edges[i][2]

            fo2 =fo2 + (((fo2_s-1)*mis[edges[i][1], edges[i][2]] )+ pis[edges[i][1], edges[i][2]])
            pis_t = pis[edges[i][1], edges[i][2]]
            mis_t = mis[edges[i][1], edges[i][2]]
            println("Corte ($a1, $a2) t: $fo2_s FO2: $fo2 , pis: $pis_t mis: $mis_t")
        end

        for i in 1:nVertex
            for j in 1:nVertex
                if i != j
                    fo3 = fo3 + (JuMP.value(X[i,j,1]) * ((vertex[i][1]-0)^2+(vertex[i][2]-0)^2 ))^(1/2)
                end
            end
        end
        #fo1 = sum(X[i, j, t] * mis[i, j] for i in 1:nVertex for j in 1:nVertex, t in 1:2*nEdges if i != j && (Ac[i,j]==0 || Ac[j,i]==0))
        #fo2 = sum(((sum( for t in 1:2*nEdges)-1)
        #            *mis[edges[i][1], edges[i][2]])+pis[edges[i][1], edges[i][2]] for i in 1:nEdges)
        #fo3 = sum(X[i,j,1] * (((vertex[i][1]-0)^2+(vertex[i][2]-0)^2 )^(1/2)) for i in 1:nVertex for j in 1:nVertex if i != j)

        println("F1: $fo1 F2: $fo2 F3: $fo3")

        popfirst!(file)
    end
end

modelo_final(lista1)
