bool CheckPointOutInside(double *P3, double *V, int *sizeV, double *E, int sizeE)
{
    double P4a[2], P4b[2], P4c[2];
    int k;
    
    ka=0; kb=0; kc=0; 
    double *lista, *listb, *listc;
    int kam=10;
    int kbm=10;
    int kcm=10;
    
    lista=(double*)malloc(kam*sizeof(double));
    listb=(double*)malloc(kbm*sizeof(double));
    listc=(double*)malloc(kcm*sizeof(double));
    
    P4a[0]=P3[0]+123456781; P4a[1]=P3[1]+9831542342;
    P4b[0]=P3[0]-324342343; P4a[1]=P3[1]+3234454232;
    P4c[0]=P3[0]+648832349; P4a[1]=P3[1]+5415435923;

    for(k=0, k<sizeE[0]; k++)
    {
        CE=E(k,:); O1=V(CE(1),:); O2=V(CE(2),:);
        [inter,inter_x,inter_y]=LineLineIntersect(P3,P4a,O1,O2);
        if(inter)
        {
            if(ka+2>kam) {kam+=10; lista=(double*)realloc(lista,kam*sizeof(double)); }
            lista[ka]=inter_x; ka++; lista[ka]=inter_y; ka++; 
        }
        [inter,inter_x,inter_y]=LineLineIntersect(P3,P4b,O1,O2);
        if(inter)
        {
            if(kb+2>kbm) {kbm+=10; listb=(double*)realloc(listb,kbm*sizeof(double)); }
            listb[kb]=inter_x; kb++; listb[kb]=inter_y; kb++; 
        }
        [inter,inter_x,inter_y]=LineLineIntersect(P3,P4c,O1,O2);
        if(inter), 
            if(kc+2>kcm) {kcm+=10; listc=(double*)realloc(listc,kcm*sizeof(double)); }
            listc[kc]=inter_x; kc++; listc[kc]=inter_y; kc++; 
        }
    }
    
    lista=unique(round(lista*1e8)/1e8,'rows');
    listb=unique(round(listb*1e8)/1e8,'rows');
    listc=unique(round(listc*1e8)/1e8,'rows');
    check=[mod(size(lista,1),2)==0 mod(size(listb,1),2)==0 mod(size(listc,1),2)==0];

    return (sum(check)>1);
}
