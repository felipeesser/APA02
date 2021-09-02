Trabalho 02 de APA

para gerar os registros das tempos de execução e gerar um caso de teste basta compilar e executar os arquivos "Floyd.c" e "TG.c":

gcc TG.c Floyd.c -lm -g -o main  
./main

Os arquivos contidos na pasta "dados_gerados" do repositório correspondem ao registro dos tempos de execução  do algoritmo de Floyd com grafos de 100 a 490 nós, tamanho de passo igual a 10 nós por medição, e com as seguintes densidades: 10%,20%,40%,50%,70%,90%.
Os gráficos foram gerados com base nos arquivos contidos na pasta "dados_gerados" e estão localizados na pasta "graficos_gerados".
Para gerar os gráficos a ferramenta gnuplot foi utilizada com os seguintes comandos:

plot "densidadeFixa(0.1)"with lines lt rgb "red","densidadeFixa(0.9)" with lines lt rgb "green","densidadeFixa(0.4)" with lines lt rgb "blue"

plot "densidadeFixa(0.1)"with lines lt rgb "red","densidadeFixa(0.9)" with lines lt rgb "brown","densidadeFixa(0.7)" with lines lt rgb "green"

plot "densidadeFixa(0.1)"with lines lt rgb "red"

plot [0:40] (x**2)

plot [0:40] (x**3)
