#include "Lista.h"
No* insereNo(No* l,int x){
    No* aux=(No*)malloc(sizeof(No));
    aux->n=x;
    aux->prox=l;
    return aux;
}
void liberalista(No *l){
    if (l)
    {   No* temp=l->prox;
        free(l);
        liberalista(temp);
    }   
}
void imprimelista(No *l){
    if (l)
    {   No* temp=l->prox;
        imprimelista(temp);
        printf("%d ",l->n);
    }   
}