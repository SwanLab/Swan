function BasisINT =  BasisIntf_projects(NDISPpr,Ui_Si,WeightedSingVal,nMODESint,f1,f2)
iini = 1 ;
BasisINT = [] ;

for ipro = 1:length(NDISPpr)
    ifin = (iini-1) + NDISPpr(ipro) ;
    BasisUdef_loc = Ui_Si(:,iini:ifin) ;
    Xa_i = [BasisUdef_loc(f1,:),BasisUdef_loc(f2,:)] ;
    [BasisINT_i S_i]= SVDT(Xa_i,0) ;
    if   WeightedSingVal == 1
        S_i = S_i/(S_i(1)) ;
        BasisINT_i = bsxfun(@times,BasisINT_i',S_i)';
    end  
    nmodesi = min(nMODESint(ipro),length(S_i));
    BasisINT = [BasisINT  BasisINT_i(:,1:nmodesi) ];
    iini = ifin+1 ;
   
end