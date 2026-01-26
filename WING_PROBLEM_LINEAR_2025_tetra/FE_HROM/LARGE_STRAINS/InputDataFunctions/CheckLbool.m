clc
clear all

CN =[1 2 5 6
    2 3 5 4] ;

nnode = 6 ;
ngaus = 2 ;
ndim = 3;


Lbool_nonvector = Lbool_nonvectorized(CN,nnode,ngaus,ndim)  ;


Lbool_vector = Lbool_vectorized(CN,nnode,ngaus,ndim) ; 


eee = full(Lbool_vector)-Lbool_nonvector