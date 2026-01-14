function MACROVAR = DispMACROtime(MACRODEF,COORrel,DATA)
%
% JAHO- 18-JAN-2O21
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')    
end


 ndim = DATA.MESH.ndim ; 
 %COOR = MESH.COOR';  
 nnode = DATA.MESH.nnode ; 
% IDENTITY_F = reshape(IDENTITY_F,ndim^2,[]);  
%  IDENTITY_F = IDENTITY_F(:,1:nnode) ; 
 % IDENTITY_F = reshape(IDENTITY_F,ndim^2,[]) ; 
%  if ndim == 2
%     IDENTITY = [1;1;0;0] ;
% elseif ndim ==3 
%     IDENTITY = [1;1;1;0;0;0;0;0;0] ; 
% end
DATA = DefaultField(DATA,'ngausTotalLOC',DATA.MESH.ngausT)  ; 
ngausTotalLOC = DATA.ngausTotalLOC ; 

nloads = length(MACRODEF) ; 
U =  zeros(nnode*ndim,nloads) ;  
a = zeros(nloads,length(DATA.STEPS)) ; 
Uf =  zeros(ngausTotalLOC*ndim^2,nloads) ;  




for  iload = 1:nloads    
    GRADuMACRO = MACRODEF(iload).AMPLITUDE ;   
  %  F_ident = Fgrad-eye(ndim) ; 
  if ~isempty(COORrel)
    dMACROmax = GRADuMACRO*COORrel ;     
    dMACROmax = dMACROmax(:) ;   
  end
    
    GRADuMACROst =  Tensor2VoigtF(GRADuMACRO)  ; 
    Uf(:,iload) = repmat(GRADuMACROst,ngausTotalLOC,1) ;   
    if ~isempty(COORrel)
    U(:,iload) =dMACROmax ;     
    end
    FactorSteps =  MACRODEF(iload).TIMEFUN(DATA.STEPS) ;
    LIMITS_interval =  MACRODEF(iload).INTERVAL ;
    FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
    a(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
end

if ~isempty(COORrel)
MACROVAR.DISP.U = U ; 
MACROVAR.DISP.a = a ; 
end
MACROVAR.GRADuMACROst.U = Uf ;
MACROVAR.GRADuMACROst.a = a ;

 


 