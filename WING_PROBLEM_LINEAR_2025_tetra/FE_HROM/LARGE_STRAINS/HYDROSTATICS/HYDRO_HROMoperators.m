function [OPERHROMhydro,DATAHROM] = HYDRO_HROMoperators(OPERFEhydro,DATAFE,DATAHROM,ECMdata_press,BasisUall)
% HROM  counterparts of hydrodynamic operators
% JAHO, 6-July-2021
if nargin == 0
    load('tmp.mat')
end

%% -------------------------------------------------------

% OPERFEhydro.tTANGiniST
% OPERFEhydro.BstB
% OPERFEhydro.Lbool
% OPERFEhydro.mNORMiniST
% OPERFEhydro.NbST
% OPERFEhydro.NbST_w
% OPERFEhydro.NbSTgrav
% OPERFEhydro.wST
% OPERFEhydro.yPRESSiniST

% DATAFE.MESH.HYDRO.nelemB ;
% DATAFE.MESH.HYDRO.ngausT ;


% ----------------------------------------------
if ~isempty(ECMdata_press)
setPoints = ECMdata_press.setPoints ;

OPERHROMhydro.wST = ECMdata_press.wRED ;
OPERHROMhydro.JacobianWeights = OPERFEhydro.JacobianWeights(setPoints) ;


OPERHROMhydro.NbSTgrav = OPERFEhydro.NbSTgrav(setPoints,:) ;
OPERHROMhydro.yPRESSiniST = OPERFEhydro.yPRESSiniST(setPoints) ;
DATAHROM.MESH.HYDRO.ngausT = length(OPERFEhydro.wST) ;
DATAHROM.MESH.HYDRO.nelemB = [] ;


ndim  =3;
INDEXloc = small2large(setPoints,ndim) ;

for idim = 1:2
    OPERHROMhydro.tTANGiniST{idim} = OPERFEhydro.tTANGiniST{idim}(INDEXloc) ;
    OPERHROMhydro.BstB{idim} = OPERFEhydro.BstB{idim}(INDEXloc,:) ;
    
end
OPERHROMhydro.mNORMiniST = OPERFEhydro.mNORMiniST(INDEXloc) ;
OPERHROMhydro.NbST = OPERFEhydro.NbST(INDEXloc,:) ;





 
nrows_perGAUSS = DATAFE.MESH.HYDRO.nnodeE*ndim ; 
INDEXloc = small2large(setPoints,nrows_perGAUSS) ;


 OPERHROMhydro.Lbool = OPERFEhydro.Lbool(INDEXloc,:)*BasisUall ;


% OPERHROMhydro.NbST_w  =  OPERFEhydro.NbST_w(INDEXloc,:)*BasisUall ;

%  This operator is to be computed again, using OPERHROMhydro.wST

% 1. Multiply each column of OPERHROMhydro.NbST by the corresponding reduced weight
%  ------------------------------------------------------------------------
NbST_w_elem =  zeros(size(OPERHROMhydro.NbST)) ; 
for idim = 1:ndim 
    NbST_w_elem(idim:ndim:end,:) = bsxfun(@times,OPERHROMhydro.NbST(idim:ndim:end,:),OPERHROMhydro.wST)  ; 
end
OPERHROMhydro.NbST_w =  ConvertBlockDiag_general(NbST_w_elem,ndim)*OPERHROMhydro.Lbool   ; 
 
 OPERHROMhydro.irowsNDIM =[] ; 
  OPERHROMhydro.icolsNDIM =[] ; 
   OPERHROMhydro.irows1 =[] ; 
  OPERHROMhydro.icols1 =[] ; 
  
else
    OPERHROMhydro = [] ; 
end
  

