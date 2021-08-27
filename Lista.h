#ifndef LISTA_H
#define LISTA_H
#include <stdlib.h>
#include <stdio.h>
typedef struct no{
        int n;
        struct no* prox;
}No;
No* insereNo(No* l,int x);
void imprimelista(No *l);
void liberalista(No *l);
#endif