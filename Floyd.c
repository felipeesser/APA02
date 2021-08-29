#include "TG.h"
#include "Lista.h"
#include <limits.h>
#include <math.h>
#define MAX INT_MAX

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
void _insere_arcos_grafo_completo(TG *g, int ordem, int distMinima, int distMaxima)
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
    while ( (densidade <= thresholdDensidade && arcosInseridos < totalArcos) ||
            (densidade >  thresholdDensidade && arcosRemovidos < maximoArcos - totalArcos)) {
        no1Sorteado = _num_aleatorio(1, ordem);
        no2Sorteado = _num_aleatorio(1, ordem);
        distancia = _num_aleatorio(distMinima, distMaxima);

        no1=no1Sorteado, no2=no2Sorteado;
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
        if (inseriuArco) {
            inseriuArco = 0;
            arcosInseridos++;
            printf("arcosIseridos: %ld\n", arcosInseridos);
        }
        if (removeuArco) {
            removeuArco = 0;
            arcosRemovidos++;
            printf("arcosRemovidos: %ld\n", arcosRemovidos);
        }
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

typedef struct elemmatrix{
        int n;
        No* k;
}EM;

void imprimematrizlista(int n,EM*** mc){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            if(i!=j){
            printf("%d ",i+1);
            imprimelista(mc[n][i][j].k);
            printf("%d ",j+1);
            printf("\n");
            }
        }
        
    }
    
}

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
    for (int i = 0; i < n; i++)
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
            printf("%d ",matrix[n][i][j].n);
        }
        
    }
    printf("\n");
    printf("\n");
    
}
int compara(int a,int b,int c,EM*** m,int k,int j,int i){
    if (b!=MAX && c!=MAX && (b+c)<a)
    {   
        m[k][i][j].k=m[k-1][i][k-1].k;
        m[k][i][j].k=insereNo(m[k][i][j].k,k);
        return b+c;
    }
    m[k][i][j].k=m[k-1][i][j].k;
    return a;
}
void floyd(int n,TG* g,int **c,EM*** m){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            m[0][i][j].n=c[i][j];
            m[0][i][j].k=NULL;
        }
        
    }
    for (int k = 1; k <= n; k++)
    {   
        for (int i = 0; i <n; i++)
        {
            for (int j = 0; j < n; j++)
            {   
                m[k][i][j].n=compara(m[k-1][i][j].n,m[k-1][i][k-1].n,m[k-1][k-1][j].n,m,k,j,i);
            }
        }
        
    }
    
    
}
int main(){
    TG *g=inicializa();
    int** matrix;
    EM*** M=criaMM(4);
    g=ins_no(g,1);
    g=ins_no(g,2);
    ins_um_sentido(g,1,2,7);
    g=ins_no(g,3);
    g=ins_no(g,4);
    ins_um_sentido(g,2,3,11);
    ins_um_sentido(g,3,4,5);
    ins_um_sentido(g,4,1,4);
    imprime(g);
    matrix=criaM(4,g);
    floyd(4,g,matrix,M);
    imprimem(matrix,4);
    imprimemm(M,4);
    printf("\n");
    imprimematrizlista(4,M);
    libera(g);
    liberam(matrix,4);
    liberamm(M,4);

    // TG *g = cria_grafo_aleatorio(1000,0.49,3,5,1);
    // imprime(g);

    // Para compilar tem que linkar a biblioteca libm
    // gcc TG.c Lista.c Floyd.c -lm -g -o main

    return 0;
}
