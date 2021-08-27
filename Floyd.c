#include "TG.h"
#include "Lista.h"
#include <limits.h>
#define MAX INT_MAX

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
    return 0;
}
