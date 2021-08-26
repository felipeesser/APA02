#include "TG.h"
#include <limits.h>
#define MAX INT_MAX
int** criaM(int n,TG* g){
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
int ***criaMM(int n){
    int*** M=(int***) malloc((n+1)*sizeof(int**));
    for(int i=0;i<=n;i++){
        M[i]=(int**)malloc(n*sizeof(int*));
        for(int j=0;j<n;j++){
            M[i][j]=(int*)malloc(n*sizeof(int));
        
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
void liberamm(int*** matrix,int n){
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
void imprimem(int** matrix,int n){
    for (int i = 0; i < n; i++)
    {   printf("\n");
        for (int j = 0; j < n; j++)
        {
            printf("%d ",matrix[i][j]);
        }
        
    }
    
}
void imprimemm(int*** matrix,int n){
    for (int k = 0; k <= n; k++)
    {printf("\n");
        for (int i = 0; i < n; i++)
    {   printf("\n");
        for (int j = 0; j < n; j++)
        {
            printf("%d ",matrix[k][i][j]);
        }
        
    }
    }
}
int compara(int a,int b){
    if (b>0 && a>b)
    {
        return b;
    }
    return a;
}
void floyd(int n,TG* g,int **c,int*** m){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            m[0][i][j]=c[i][j];
        }
        
    }
    for (int k = 1; k <= n; k++)
    {
        for (int i = 0; i <n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                m[k][i][j]=compara(m[k-1][i][j],(m[k-1][i][k]+m[k-1][k][j]));
            }
            
        }
        
    }
    
    
}
int main(){
    TG *g=inicializa();
    int** matrix;
    int*** M=criaMM(4);
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
    libera(g);
    liberam(matrix,4);
    liberamm(M,4);
    return 0;
}
