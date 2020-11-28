using LinearAlgebra
function lsqm(vx,vy) #Regressão linear usada para os dois exercícios
    m = size(vx',2)
    sumyi, sumxi, sumxiyi, sumxisq, sumyisq = 0, 0, 0, 0, 0
    for i = 1:m
        sumyi = sumyi + vy'[i]
        sumxi = sumxi + vx'[i]
        sumxiyi = sumxiyi + vx'[i] * vy'[i]
        sumxisq = sumxisq + vx'[i] ^ 2
        sumyisq = sumyisq + vy'[i] ^ 2
    end
    α = m * sumxisq - sumxi ^ 2  
    if α != 0
        a₀ = 1/α * (sumyi * sumxisq - sumxi * sumxiyi)
        a₁ = 1/α * (m * sumxiyi - sumxi * sumyi)
    end
    
    #=for j = 2011:2090
        f(j) = a₀ + a₁*j
        k = f(j) - 13.72
        println("Temperatura média global no ano ", j ," é de: ",f(j))
        if k < 2.0
            println("A meta estabelecida pelo Acordo de Paris foi cumprida: mθ(x) - 13,72 = ", k)
        elseif k ≥ 2.0
            println("A meta estabelecida pelo Acordo de Paris não foi cumprida: mθ(x) - 13,72 = ", k)
        end
    end =#  
end

#=== Exercício 1 ===#
#===   Letra a)  ===#
#= Para obter a % de Pessoas com Ensino Superior no estado i, foi construída a seguinte matriz A, com dados de cidades.ibge.gov.br, na seção Pesquisas > Censo > 
  Amostra-Religião > Pessoas com 25 anos ou mais de idade > Nível de instrução, em que as linhas de A são os estados, na ordem da tabela e cada coluna de A é um nível de instrução diferente, na ordem do site do IBGE=#
  #= A = [ 189839      41288      76511      30434      813 
           1022769     180117     290571     110974     2758
           128900      41881      101143     33198      997
           778066      223462     479657     133539     8714
           4540315     911715     1824996    499196     18297
           2624110     618656     1013821    329085     8653
           1026023     298637     520078     229988     3862
           1734991     525522     848600     356609     8291
           1967557     382793     669941     173766     6868
           871746      252095     380247     176549     5209
           704713      199879     318539     166902     1680
           6298868     1637102    2547479    1242682    25239
           2083485     547090     822764     228974     7177
           1317438     222898     400653     169584     3085
           3059025     936133     1423741    793391     12000
           2795033     600645     1108466    392863     10616
           1079572     190157     289370     122788     2742
           3807398     1735282    3023446    1435229    28726
           995748      217127     412887     147701     2512
           3215968     1110416    1622785    758000     12051
           476834      112314     174529     66856      1341
           96418       27620      64999      21438      522
           1744626     625598     924654     473102     6334
           10437806    4135385    6932185    3843068    109536
           637706      133442     241307     94616      1954
           380145      86801      172848     73210      1307] 
     e usada a seguinte função:
     function ex1(A)
         k = size(A,1)
      y, z = zeros(k), zeros(k)
      for i = 1:size(A,1)
            for j = 1:size(A,2)
               y[i] += A[i,j]
            end
      end 
      z[:] = A[:,4]
      z[:] = z[:]./y[:]
      #a = round.(z[:].*100, digits = 2)
      return z # retorna a % de Pessoas com Ensino Superior nos estados
end =#
# Para calcular o PIB per capita:
#= Pibrs = [ 8477
          24575
          8266
          59779
          154340
          77865
          82122
          97576
          45256
          59600
          43514
          351381
          77848
          31947
          217290
          95187
          22060  
          407123
          32339
          252483
          23561
          6341
          152482
          1247596
          23932
          17240 ]; # Vetor que contém o pib em R$ dos estados na ordem da tabela.
 P = [733559
      3120494
      669526
      3483985
      14016906
      8452381 
      3514952
      6003788
      6574789
      3035122
      2449024
      19597330
      7581051
      3766528
      10444526
      8796448
      3118360
      15989929
      3168027
      10693929
      1562409
      450479
      6248436
      41262199
      2068017
      1383445 ] # Vetor que contém a população dos estados na ordem da tabela =#

      #PC = (Pibrs[:]./P[:])
      #return PC #retorna o PIB per capita em R$ (não em mil R$ como no site que foi retirado o vetor Pibrs)
      #= PC arredondado:  
       PC= [0.012
            0.008
            0.012
            0.017
            0.011
            0.009
            0.023
            0.016
            0.007
            0.02
            0.018
            0.018
            0.01
            0.008
            0.021
            0.011
            0.007
            0.025
            0.01
            0.024
            0.015
            0.014
            0.024
            0.03
            0.012
            0.012]=#
      #Observação: para o calculo do modelo foi utilizado os valores não arredondados de PC 


#===   Letra b)  ===#
#=scatter((z[:],PC[:]), label="Dados", c=:skyblue, xlims=(4,16), ylims=(500,3200), xlabel ="% de Pessoas com ensino superior no estado i", ylabel = "PIB per capita do Estado i", legend=:topleft)=#

#===   Letra c)  ===#
# Fazendo θ₀, θ₁ = lsqm(z, PC) temos que:
# mθ(x) = θ₀ + θ₁*x é dado por (θ₀, θ₁) = (-653.6509770358325, 227.01163531896373)
#=scatter((z[:],PC[:]), label="Dados", c=:skyblue, xlims=(4,16), ylims=(500,3200), xlabel ="% de Pessoas com ensino superior no estado i", ylabel = "PIB per capita do Estado i", legend=:topleft)=#
#plot!(x-> a0 + a1*x, 0, 20, label="Projeção Linear", c=:red)

#===   Letra d)  ===#
# Percebe-se, pelo gráfico, que quanto maior a porcentagem de pessoas com ensino superior em um estado, maior seu PIB per capita, assim como, quanto menor for essa porcentagem, menor o PIB per capita.

#=== Exercício 2 ===#
#===   Letra a)  ===#
# Dados:
#= vx=[1970
   1971
   1972
   1973
   1974
   1975
   1976
   1977
   1978
   1979
   1980
   1981
   1982
   1983
   1984
   1985
   1986
   1987
   1988
   1989
   1990
   1991
   1992
   1993
   1994
   1995
   1996
   1997
   1998
   1999
   2000
   2001
   2002
   2003
   2004
   2005
   2006
   2007
   2008
   2009
   2010]

 vy=[14.030
   13.900
   14.000
   14.140
   13.920
   13.950
   13.840
   14.120
   14.010
   14.080
   14.190
   14.260
   14.040
   14.250
   14.090
   14.040
   14.120
   14.270
   14.310
   14.190
   14.360
   14.350
   14.130
   14.130
   14.230
   14.370
   14.290
   14.390
   14.560
   14.320
   14.330
   14.470
   14.560
   14.550
   14.480
   14.620
   14.550
   14.580
   14.440
   14.580
   14.630] =#
#scatter(vx[:],vy[:], color=:skyblue, label=false)
# θ₀, θ₁ = lsqm(vx,vy)
# scatter((vx[:],vy[:]), label="Dados", c=:skyblue, xlabel = "Ano", ylabel="Temperatura média global", xlims=(1965,2015), ylims=(13.5,15))
# plot!(x-> θ₀ + θ₁*x, 0, 2500, color=:red, label=Projeção Linear, legend=:topleft)
#= Veja que para melhorar a visualização da regressão linear feita, foi aplicada a visualização somente no intervalo x = [1965, 2015] e y = [13,5, 15],
e o modelo de regressão linear mθ(x) = θ₀ + θ₁*x é dado por (θ₀, θ₁) = (-18,897177700296563, 0,01666202090590055). =#

#===   Letra b)  ===#
# Seja mθ(x) = -18,897177700296563 + 0,01666202090590055*x, queremos calcular mθ(2090), cujo resultado é mθ(2090) = 15,926 graus Celsius.
# mθ(2090) - 13,72 = 2,206, que é superior a meta estabelecida no  Acordo de Paris.

#===   Letra c)  ===#
# A projeção feita no item b sugere que a temperatura média global é crescente (de acordo com a projeção com base em dados anteriores).
# Veja que, de 2011 a 2090 podemos calcular a projeção de acordo com o modelo feito e verificar em qual ano o Acordo de Paris não se cumpre e que a temperatura média global é de fato crescente.
# Passo 1:  i = 2011, ..., 2090
# Passo 2:  Calcule mθ(i)
# Passo 3:  se mθ(i) - 13,72 < 2.0 então a meta foi cumprida
# Passo 4:  se mθ(i) - 13,72 ≥ 2.0 então a meta não foi cumprida, retorne i
# Dessa forma, verificamos que a partir do ano 2078 a meta não será cumprida, com mθ(x) - 13,72 = 2,006.
