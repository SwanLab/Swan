function [G,uBAR,DOFr,DOFm,AREA,R] = ...
    PRESCRIBED_FLUCTUATIONS_J1face_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA,FIXED_FACE)

if nargin == 0
    load('tmp1.mat')
end

% Loading information of fluctuations (of one of the faces)
load(DATA.FLUCTUATIONS_BOUNDARY_NameWS,'RIGID_BODY_MATRIX','FLUCTUATION_MATRIX','MASS_MATRIX_GEOMETRIC')  ;
% -----------------------------------------------------
ndim = 3;
AREA = zeros(2,1) ;
%%%% FACE 1
iface=1 ;
R= [] ; 
[DOFrA,DOFmA,JA,uA,AREAA,RA,MstA] = VariablesFaces(a_A,iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,FLUCTUATION_MATRIX,RIGID_BODY_MATRIX,...
    MASS_MATRIX_GEOMETRIC,FIXED_FACE) ;



%%%% FACE 2
iface=2 ;
[DOFrB,DOFmB,JB,uB,AREAB,RB,MstB] = VariablesFaces(a_B,iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,FLUCTUATION_MATRIX,RIGID_BODY_MATRIX,...
    MASS_MATRIX_GEOMETRIC,FIXED_FACE) ;

AREA = [AREAA,AREAB]; 
R = {RA,RB}; 



DOFr = [DOFrA; DOFrB] ;
DOFm = [DOFmA; DOFmB] ;


G = sparse(length(DOFr),length(DOFm)) ;
uBAR = zeros(length(DOFr),1) ;


% COLUMN 1 
% --------- 
if isempty(DOFmA)
    
     iini = 1;
    ifin = length(DOFrA) ;  
     uBAR(iini:ifin) = uA  ;
    
    iini = ifin +1;
    ifin = length(DOFr) ;  
    
    G(iini:ifin,:) = JB ;
    uBAR(iini:ifin) = uB  ;
else
     iini = 1;
    ifin = length(DOFrA) ;  
    
    G(iini:ifin,:) = JA ;
    uBAR(iini:ifin) = uA  ;
    
        iini = ifin +1;
    ifin = length(DOFr) ;  
     uBAR(iini:ifin) = uB  ;
    
end
  


end 

function  [DOFrA,DOFmA,JA,uA,AREA,R,Mst] = VariablesFaces(a_A,iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,FLUCTUATION_MATRIX,RIGID_BODY_MATRIX,...
    MASS_MATRIX_GEOMETRIC,FIXED_FACE)

RIGID_BODY_AMPLITUDE = a_A ;
if FIXED_FACE == iface
     [DOFrA,DOFmA,JA,uA,AREA,R,Mst] ...
    =  FLuctuBC_oneface_FLZERO(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    RIGID_BODY_AMPLITUDE) ;
else
    BasisINTfluct = FLUCTUATION_MATRIX{1} ;
    BasisINTrb = RIGID_BODY_MATRIX{1} ;
    M = MASS_MATRIX_GEOMETRIC{1} ;
    if iscell(M)
        M = M{1} ;
        BasisINTrb = BasisINTrb{1} ;
    end
    [DOFrA,DOFmA,JA,uA,AREA,R,Mst] ...
    =  FLuctuBC_oneface(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,BasisINTfluct,...
    RIGID_BODY_AMPLITUDE,BasisINTrb,M) ;
end


end

