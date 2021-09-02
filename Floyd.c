#include "TG.h"
#include <limits.h>
#include <time.h> 
#include <math.h>
#include <sys/resource.h>
#define MAX INT_MAX

void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
{
  long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
  struct rusage ptempo;

  getrusage(0,&ptempo);

  seg_CPU = ptempo.ru_utime.tv_sec;
  mseg_CPU = ptempo.ru_utime.tv_usec;
  seg_sistema = ptempo.ru_stime.tv_sec;
  mseg_sistema = ptempo.ru_stime.tv_usec;

 *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
 *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
}
/**
Mapeia uma contagem linear em uma contagem circular; correspondendo o menor
elemento da contagem linear ao menor elemento da contagem circular.

inicioL = menor valor da sequencia linear
inicioC = menor valor da sequência circular
fimC = maior valor da sequência circular

ex.:
s = _map_circular(e,inicioL=0,inicioC=1,fimC=3);

-------
e  |  s
-------
0  |  1
1  |  2
2  |  3
3  |  1
4  |  2
5  |  3
...|...
-------
(e)ntrada|(s)aida
*/
int _map_circular(unsigned int entradaL, unsigned int inicioL, unsigned int inicioC, unsigned int fimC)
{
    return (entradaL-inicioL)%(fimC-inicioC+1)+inicioC;
}

int _num_aleatorio(int valorMinimo, int valorMaximo)
{
    return (rand() % (valorMaximo - valorMinimo + 1)) + valorMinimo;
}

double _clamp(double entrada, int valorMinimo, int valorMaximo)
{
    double aux = entrada < valorMinimo ? valorMinimo : entrada;
    return aux > valorMaximo ? valorMaximo : aux;
}

/**
Recebe um grafo e insere arcos com a finalidade de transformar o grafo em um
grafo completo;
*/
void _insere_arcos_grafo_completo(TG *g, size_t ordem, int distMinima, int distMaxima)
{
    int distancia=0;
    size_t i=0, j=0;
    for (i=0; i<ordem; i++) {
        for (j=0; j<ordem; j++) {
            distancia = _num_aleatorio(distMinima, distMaxima);
            if (i != j) {
                ins_arco(g, i+1, j+1, distancia);
            }
        }
    }
}

/**
Insere ou remove arcos para que o grafo tenha a densidade desejada;

Essa função deve ser usada em conjunto com a função _insere_arcos_grafo_completo
*/
void _insere_remove_arcos(TG *g, size_t ordem, double densidade, double thresholdDensidade, int distMinima, int distMaxima)
{
    size_t totalArcos=(size_t) ceil(ordem*(ordem-1)*densidade);
    size_t maximoArcos=ordem*(ordem-1);
    int distancia=0;
    int no1Sorteado=0, no2Sorteado=0, no1=0, no2=0;
    int inseriuArco=0, removeuArco=0;
    size_t i=0, j=0, arcosInseridos=0, arcosRemovidos=0;

    // NOTE - itera enquanto ainda for preciso inserir ou remover nos
    while ( (densidade <= thresholdDensidade && arcosInseridos < totalArcos) ||
            (densidade >  thresholdDensidade && arcosRemovidos < maximoArcos - totalArcos)) {
        // NOTE - sorteia um arco a ser inserido ou alterado
        no1Sorteado = _num_aleatorio(1, ordem);
        no2Sorteado = _num_aleatorio(1, ordem);
        distancia = _num_aleatorio(distMinima, distMaxima);

        no1=no1Sorteado, no2=no2Sorteado;
        // NOTE - os dois fors são para iterar por todos os possiveis arcos do
        // grafo caso o arco sorteado já tenha sido inserido ou removido
        // NOTE - iteração dos fors continua até que um arco seja inserido ou
        // removido
        for (i=0; i<ordem && (!inseriuArco && !removeuArco); i++) {
            no1 = _map_circular(no1Sorteado+i, 1, 1, ordem);
            for (j=0; j<ordem && (!inseriuArco && !removeuArco); j++) {
                no2 = _map_circular(no2Sorteado+j, 1, 1, ordem);
                if (no1 != no2) {
                    if (densidade > thresholdDensidade) {
                        removeuArco = retira_arco(g, no1, no2);
                    } else {
                        inseriuArco = ins_arco(g, no1, no2, distancia);
                    }
                }
            }
        }
        // NOTE - caso tenha inserido um arco, aumenta o contador de arcos inseridos
        // NOTE - caso tenha inserido um arco, aumenta o contador de arcos inseridos
        if (inseriuArco) {
            inseriuArco = 0;
            arcosInseridos++;
            // printf("arcosIseridos: %ld\n", arcosInseridos);
        }
        else if (removeuArco) {
            removeuArco = 0;
            arcosRemovidos++;
            // printf("arcosRemovidos: %ld\n", arcosRemovidos);
        }
        // NOTE - volta e sorteia mas um arco para ser inserido
    }
}

/**
Retorna um ponteiro para um grafo criado aleatóriamente de ordem, densidade e
peso dos arcos entre distMinima e distMaxima;
*/
TG *cria_grafo_aleatorio(size_t ordem, double densidade, int distMinima, int distMaxima, int seed)
{
    densidade = _clamp(densidade, 0.0, 1.0);
    TG *novo = inicializa();
    const double thresholdDensidade = 0.6;
    int distancia=0;
    size_t i=0, j=0;
    srand(seed);

    // NOTE - adiciona os nós ao grafo
    for (i=0; i<ordem; i++) {
        // rotulos começam em 1
        novo = ins_no(novo, i+1);
    }

    // NOTE - faz um grafo completo se a densidade for maior que o threshold
    if (densidade > thresholdDensidade) {
        _insere_arcos_grafo_completo(novo, ordem, distMinima, distMaxima);
    }

    // NOTE - de acordo com o threshold, insere arcos ou remove arcos até atender a densidade
    _insere_remove_arcos(novo, ordem, densidade, thresholdDensidade, distMinima, distMaxima);

    return novo;
}

/**
Tupla para armazenar a distância e o nó usada para o caminho entre dois nós;
*/
typedef struct elemmatrix{
        int distancia;
        size_t no;
}EM;

int** criaM(int n,TG* g){//matriz inicial do grafo
    int** matrix=(int**) malloc(n*sizeof(int*));
    for(int i=0;i<n;i++){
        matrix[i]=(int*)malloc(n*sizeof(int));
    }
    TG*aux= g;
    int achou=0;
    TVIZ*aux2;
    while (aux){ 
        for (int j = 1; j <= n; j++)
        {   
            if(j==aux->id_no){
                matrix[j-1][j-1]=0;
            }
            else{
                aux2=aux->prim_viz;
                achou=0; 
                while (aux2 && !achou)
                {  
                    if(j==aux2->id_viz){
                        matrix[aux->id_no-1][j-1]=aux2->distancia;
                        achou=1;
                    
                    }
                    aux2=aux2->prox_viz;
                }

                if(!achou){
                    matrix[aux->id_no-1][j-1]=MAX;
                }
            }  
        }
        aux=aux->prox_no;
    }
    return matrix;
}
EM ***criaMM(int n){
    EM*** M=(EM***) malloc((n+1)*sizeof(EM**));
    for(int i=0;i<=n;i++){
        M[i]=(EM**)malloc(n*sizeof(EM*));
        for(int j=0;j<n;j++){
            M[i][j]=(EM*)malloc(n*sizeof(EM));
        
        }
    }
    return M;
}
void liberam(int** matrix,int n){
    for (int i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}
void liberamm(EM*** matrix,int n){
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}
void imprimem(int** matrix,int n){//k=0
    for (int i = 0; i < n; i++)
    {   printf("\n");
        for (int j = 0; j < n; j++)
        {
            printf("%d ",matrix[i][j]);
        }
        
    }
    printf("\n");
    printf("\n");
    
}
void imprimemm(EM*** matrix,int n){//k=n
        for (int i = 0; i < n; i++)
    {   printf("\n");
        for (int j = 0; j < n; j++)
        {
            int dist = matrix[n][i][j].distancia;
            if (dist != INT_MAX) {
                printf("%2d ", dist);
            } else {
                printf("oo ");
            }
        }
    }
    printf("\n");
    printf("\n");
    
}

/**
Retorna o mínimo entre {ij, ik+kj}
*/


int compara(int ij, int ik, int kj)
{
    if (ik == INT_MAX || kj == INT_MAX) {
        return ij;
    }
    else if (ik > 0 && kj > INT_MAX - ik) {
        /* overflow */
        printf("int não comporta a soma das distancias\n");
        exit(1);
    } else if (ik < 0 && kj < INT_MIN - ik) {
        /* underflow */
        printf("int não comporta a soma das distancias\n");
        exit(1);
    }
    return ij < ik+kj ? ij : ik+kj;
}


/**
Função interna ao imprime_caminho que percorre recursivamente os caminhos, e
imprime cada passo.
*/
void _passa_por(EM ***c, size_t ordem, size_t i, size_t j)
{
    int passaPor = c[ordem][i][j].no;
    if (passaPor==0) {
        return;
    } else if (i == passaPor-1 || j == passaPor-1) {
        printf("(entra em loop...) - ");
        return;
    } else {
        _passa_por(c, ordem, i, passaPor-1);
        printf("%d - ", passaPor);
        _passa_por(c, ordem, passaPor-1, j);
    }
}
/**
Imprime todos os caminhos de menor custo encontrados.

Recebe a matriz gerada pelo algoritmo de Floyd
*/
void imprime_caminho(EM ***resFloyd, size_t ordem)
{
    size_t i=0,j=0;
    for (i=0; i<ordem; i++) {
        for (j=0; j<ordem; j++) {
            if (resFloyd[ordem][i][j].distancia != INT_MAX && resFloyd[ordem][i][j].distancia != 0) {
                printf("De %ld até %ld (custo: %3d): %ld - ",i+1 ,j+1, resFloyd[ordem][i][j].distancia, i+1);
                _passa_por(resFloyd,ordem,i,j);
                printf("%ld\n", j+1);
            }
        }
    }
}


void floyd(int n,TG* g,int **c,EM*** m){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            m[0][i][j].distancia=c[i][j];
            m[0][i][j].no=0;
        }
        
    }
    for (int k = 1; k <= n; k++)
    {
        for (int i = 0; i <n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int dist_ij = m[k-1][i][j].distancia;
                m[k][i][j].distancia=compara(dist_ij, m[k-1][i][k-1].distancia, m[k-1][k-1][j].distancia);
                if (m[k][i][j].distancia != dist_ij) {
                    m[k][i][j].no = k;
                } else {
                    m[k][i][j].no = m[k-1][i][j].no;
                }
            }
        }
    }
}




int main(){
// TESTE ---------------------
    TG *g;
    int** matrisI;
    EM*** M;
    g=inicializa();
    M=criaMM(4);
    g=ins_no(g,1);
    g=ins_no(g,2);
    ins_um_sentido(g,1,2,7);
    g=ins_no(g,3);
    g=ins_no(g,4);
    ins_um_sentido(g,1,3,30);
    ins_um_sentido(g,2,3,11);
    ins_um_sentido(g,3,4,5);
    ins_um_sentido(g,4,1,4);
    imprime(g);
    matrisI=criaM(4,g);
    floyd(4,g,matrisI,M);
    imprimem(matrisI,4);
    imprimemm(M,4);
    printf("\n");
    imprime_caminho(M,4);
    libera(g);
    liberam(matrisI,4);
    liberamm(M,4);
// ---------------------------

// // grafo do slide 217 --------
//     const int ordem = 3;
//     TG *g = inicializa();
//     g = ins_no(g,1);
//     g = ins_no(g,2);
//     g = ins_no(g,3);
//     ins_arco(g,1,2,4);
//     ins_arco(g,2,1,6);
//     ins_arco(g,1,3,11);
//     ins_arco(g,3,1,3);
//     ins_arco(g,2,3,2);

//     int **A0 = criaM(ordem,g);
//     EM ***A = criaMM(ordem);

//     floyd(3,g,A0,A);
//     imprimem(A0,ordem);
//     imprimemm(A,ordem);
//     imprime_caminho(A,ordem);
//     libera(g);
//     liberam(A0, ordem);
//     liberamm(A, ordem);
// // ---------------------------

// // grafo qualquer ------------
//     const int ordem = 4;
//     TG *g = inicializa();
//     g = ins_no(g,1);
//     g = ins_no(g,2);
//     g = ins_no(g,3);
//     g = ins_no(g,4);
//     ins_arco(g,1,2,2);
//     ins_arco(g,1,4,2);
//     ins_arco(g,2,3,10);
//     ins_arco(g,2,4,2);
//     ins_arco(g,4,3,2);

//     int **A0 = criaM(ordem,g);
//     EM ***A = criaMM(ordem);

//     floyd(ordem,g,A0,A);
//     imprimem(A0,ordem);
//     imprimemm(A,ordem);
//     imprime_caminho(A,ordem);
//     libera(g);
//     liberam(A0, ordem);
//     liberamm(A, ordem);
// // ---------------------------

// // grafo qualquer (ciclo negativo) ------------
//     const int ordem = 5;
//     TG *g = inicializa();
//     g = ins_no(g,1);
//     g = ins_no(g,2);
//     g = ins_no(g,3);
//     g = ins_no(g,4);
//     g = ins_no(g,5);
//     ins_arco(g,1,2,1);
//     ins_arco(g,1,3,1);
//     ins_arco(g,3,2,1);
//     ins_arco(g,2,4,4);
//     ins_arco(g,4,3,-6);
//     ins_arco(g,4,5,1);

//     int **A0 = criaM(ordem,g);
//     EM ***A = criaMM(ordem);

//     floyd(ordem,g,A0,A);
//     imprimem(A0,ordem);
//     imprimemm(A,ordem);
//     imprime_caminho(A,ordem);
//     libera(g);
//     liberam(A0, ordem);
//     liberamm(A, ordem);
// // ---------------------------

    // EM ***m = criaMM(3);
    // liberamm(m, 3);

    // TG *g = cria_grafo_aleatorio(1000,0.49,3,5,1);
    // imprime(g);

    double s_CPU_inicial,s_CPU_final,s_total_inicial,s_total_final;
    double tempo;
    FILE *f;
    size_t iTamanho, passoTamanho = 10, maxTamanho = 500;
    double iDensidade = 0.025, passoDensidade = 0.025, maxDensidade = 1.0;
    size_t iRepeticao;
    // NOTE - para cada tamanho {100, ..., }
    //double fixadensidade=0.2;
    //double fixadensidade=0.5;
    double fixadensidade=0.7;
    //double fixadensidade=0.9;
    //double fixadensidade=0.1;
    //double fixadensidade=0.4;

    for (iTamanho = 100; iTamanho < maxTamanho; iTamanho += passoTamanho)
    {   printf("%ld\n",iTamanho);
        // densidade fixa
            // NOTE - repita 10 vezes
            tempo=0;
            s_CPU_inicial=0,s_CPU_final=0,s_total_inicial=0,s_total_final=0;
            for (iRepeticao = 0; iRepeticao < 10; iRepeticao++)
            {   
                g=inicializa();
                g=cria_grafo_aleatorio(iTamanho,fixadensidade,1,100,5);
                matrisI=criaM(iTamanho,g);
                M=criaMM(iTamanho);
                Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
                floyd(iTamanho,g,matrisI,M);
                Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
                tempo+=s_CPU_final - s_CPU_inicial;
                libera(g);
                liberam(matrisI,iTamanho);
                liberamm(M,iTamanho);
            }
            f = fopen("densidadeFixa(0.7)", "a");
            fprintf (f,"%ld  %f\n",iTamanho,tempo/10);
            fclose(f);
            // TODO - calcular a media das 10 repetições
            // TODO - acrescentar a media calculada a um arquivo
    }
    size_t tamanhofixo=200;
    // NOTE - tamanho fixo
    // NOTE - para cada desisidade {10, ..., }
    for (iDensidade = 0.025; iDensidade < maxDensidade; iDensidade += passoDensidade)
    {   printf("%f\n",iDensidade);
        // NOTE - repita 10 vezes
        tempo=0;
        s_CPU_inicial=0,s_CPU_final=0,s_total_inicial=0,s_total_final=0;
        for (iRepeticao = 0; iRepeticao < 10; iRepeticao++)
        {
            g=inicializa();
            g=cria_grafo_aleatorio(tamanhofixo,iDensidade,1,100,5);
            matrisI=criaM(tamanhofixo,g);
            M=criaMM(tamanhofixo);
            Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
            floyd(tamanhofixo,g,matrisI,M);
            Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
            tempo+=s_CPU_final - s_CPU_inicial;
            libera(g);
            liberam(matrisI,tamanhofixo);
            liberamm(M,tamanhofixo);
        }
        f = fopen("tamanhoFixo(200)", "a");
        fprintf (f,"%f  %f\n",iDensidade,tempo/10);
        fclose(f);
        // TODO - calcular a media das 10 repetições
        // TODO - acrescentar a media calculada a um arquivo
    }

    return 0;
    // NOTE - para compilar o programa tem que linkar a biblioteca libm
    // gcc TG.c Lista.c Floyd.c -lm -g -o main
}
