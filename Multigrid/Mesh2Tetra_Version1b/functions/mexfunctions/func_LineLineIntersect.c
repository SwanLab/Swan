double nnz(double x)
{
    if(absd(x)<1e-16) { return 1e-16; } else { return x; }
}

void LineLineIntersect(double *P1,double *P2,double *O1,double *O2,double *Val)
{
    double inter_x, inter_y, SlopeA, SlopeB, t, disp, diso, OffsetA, OffsetB, temp;
    int swap;
    bool inter;
    bool check;
    
    check=(mind(P1[0],P2[0])>maxd(O1[0],O2[0]))||(mind(P1[1],P2[1])>maxd(O1[1],O2[1]))||(maxd(P1[0],P2[0])<mind(O1[0],O2[0]))||(maxd(P1[1],P2[1])<mind(O1[1],O2[1]));
    if(check) { Val[0]=0; Val[0]=0; Val[0]=0; return; }

    swap=0;
    
    /* Distance in x direction */
    disp=absd(P2[0]-P1[0]); diso=absd(O2[0]-O1[0]);

    /* Check if the line is not vertical, in that case swap x and y */
    if((disp<1e-10)||(diso<1e-10))
    {
        swap=1;
        t=P1[0]; P1[0]=P1[1]; P1[1]=t; t=P2[0]; P2[0]=P2[1]; P2[1]=t;
        t=O1[0]; O1[0]=O1[1]; O1[1]=t; t=O2[0]; O2[0]=O2[1]; O2[1]=t;
    }

    /* Swap to get x direction positive */
    if(P1[0]>P2[0])
    { 
        t=P1[0]; P1[0]=P2[0]; P2[0]=t; 
        t=P1[1]; P1[1]=P2[1]; P2[1]=t;  
    }
    if(O1[0]>O2[0])
    { 
        t=O1[0]; O1[0]=O2[0]; O2[0]=t; 
        t=O1[1]; O1[1]=O2[1]; O2[1]=t;  
    }
    
    /* Describe line by linear equations */
    SlopeA=(P2[1]-P1[1])/nnz(P2[0]-P1[0]); OffsetA=P1[1]-(SlopeA*P1[0]);
    SlopeB=(O2[1]-O1[1])/nnz(O2[0]-O1[0]); OffsetB=O1[1]-(SlopeB*O1[0]);

    /* Calculate intersection point */
    inter_x=(OffsetB-OffsetA) / nnz(SlopeA-SlopeB);

    SlopeA=1/nnz(SlopeA); OffsetA=P1[0]-(SlopeA*P1[1]);
    SlopeB=1/nnz(SlopeB); OffsetB=O1[0]-(SlopeB*O1[1]);
    inter_y=(OffsetB-OffsetA) / nnz(SlopeA-SlopeB);
    if(swap==1)
    {
        temp=inter_x; inter_x=inter_y; inter_y=temp;
    }

    /* Check if intersection point is between the points describing the lines */
    
    inter=(inter_x>(P1[0]-1e-15))&&(inter_x<(P2[0]+1e-15))&&(inter_x>(O1[0]-1e-15))&&(inter_x<(O2[0]+1e-15));
    
    if(inter) { Val[0]=1; } else { Val[0]=0; }
    Val[1]=inter_x;
    Val[2]=inter_y;
}
