function BND_Dyn = AccelerationDirichletBoundary(DynDATA,OPERfe,DATA)
% Dynamic problems. Acceleration 
 nDOFTOT = OPERfe.ndof ; 



if  ~isempty(DynDATA)  && DynDATA.DYNAMIC ==1  %%
    error('Option not implemented yet')
    gBOUND = DATABC(itraj).gBOUND;
    % Acceleration of BCs
   
    [d0,v0,gBOUNDdd,gBOUNDd]=AccelDirichBC(DynDATA,nDOFTOT,DOFm,DOFs,DOFl,Gbound,gBOUND,VECTOR_TIME);
else
    d0 = zeros(nDOFTOT,1) ; v0=zeros(nDOFTOT,1) ; gBOUNDdd = [] ; gBOUNDd = [] ;
end

BND_Dyn.d0  = d0 ; 
BND_Dyn.v0  = v0 ; 
BND_Dyn.gBOUNDdd  = gBOUNDdd ; 
BND_Dyn.gBOUNDd  = gBOUNDd ; 
