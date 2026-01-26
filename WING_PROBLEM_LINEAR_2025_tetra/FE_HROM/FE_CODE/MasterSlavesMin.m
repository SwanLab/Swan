function [DOFr,DOFm,G,dR,dispMACRO] = ...
    MasterSlavesMin(NODESpl,NODESpnt,CNb,COOR,TypeElementB,strainINP,DATA)
% See ALEX_THESIS_mine.pdf
%dbstop('5')
if nargin==0
    load('tmp1.mat')
end


 
% Computation of matrix R  (nstrain x nelemBND*ndime)
NODESbnd = [] ;
for i = 1:size(NODESpl,1)
    for j = 1:size(NODESpl,2)
        NODESbnd = [NODESbnd;NODESpl{i,j}] ;
    end
end
NODESbnd = unique(NODESbnd) ;


% Matrix R, such that R*dB = V*eMACRO
R = MinCond_Rmatrix(CNb,NODESbnd,COOR,TypeElementB) ;

% Master/slave nodes. G matrix

 

%

%% MASTER/slaves DOFs  %
%dbstop('32')
ndim = size(COOR,2) ;
DOFb = Nod2DOF(NODESbnd,ndim) ; % Boundary DOFs
DOFpnt = Nod2DOF(NODESpnt,ndim) ; % Corner DOFs 

% Select DOFs for eliminating rigid body motions 
 
PointsRef = NODESpnt(1:ndim) ; 
DOFp = RigidBodyPoints(PointsRef,COOR) ; 
 


DOFpnt = setdiff(DOFpnt,DOFp) ; 
[~,idx]=licols(R(:,DOFpnt)) ; 
DOFs = DOFpnt(idx) ; % Master Fixed 
DOFr = [DOFs; DOFp] ;            % Slave DOFs
DOFm = setdiff(DOFb,DOFr) ;      % Master DOFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix G 
Gs = -inv(R(:,DOFs))*R(:,DOFm) ; 
Gp = sparse(length(DOFp),size(Gs,2)) ;  
G = [Gs;Gp] ; 

% Prescribed displacement 
% -----------------------
switch DATA.REFERENCE_POINT
    case 'CORNER'
        C = COOR(NODESpnt(1),:)';
        
    case 'CENTER'  % % % 3) Displacement are measured taking as reference NODESpnt(1)
        C = (DATA.CENTER)' ;
    otherwise
        error('OPTION NOT IMPLEMENTED')
end

 
COORref =(bsxfun(@plus,COOR',-C));
%%%

% uBARp
% -----
dispMACRO = strainINP*COORref ;   % Macroscopic displacements of all DOFs
dispMACRO = dispMACRO(:); 
uBARp = dispMACRO(DOFp) ; 
% uBARs 
eMACRO = Mat2VoigtSTRAIN(strainINP) ; % Voigt notation
V = DATA.VOL_RVE ; % VOLUME 
uBARs =  -inv(R(:,DOFs))*(eMACRO*V - R(:,DOFp)*uBARp) ; 
dR = [uBARs ; uBARp] ;





 

    %





% 
% % Assembly of matrix Gi
% %----------------------
% [Gi nods  nodm DIFFcoor]= GmatrixSlaveMaster(MASTER,SLAVES,COORref')  ;
% 
% %Assembly of G
% G= sparse(3*size(Gi,1),3*size(Gi,2)) ;
% n = size(Gi,1) ; m=size(Gi,2) ;
% for i=1:3
%     ROWS = i:3:3*n ; 
%     COLS = i:3:3*m ; 
%     G(ROWS,COLS ) = Gi ;
% end
% 
% % Prescribed displacement
% % % 3) Displacement are measured taking as reference NODESpnt(1)
% % %  Accordingly, we redefine the coordinate matrix
% % 2) Center of the unit cell
% 
% 
% % Thus
% uB = MACRODEF*DIFFcoor';
% 
% % dispMACRO
% %dbstop('99')
% dispMACRO = MACRODEF*COORref ;
% dispMACRO = dispMACRO(:); 
% 
% 
% 
% DOFr = zeros(3*length(nods),1) ;
% dR = zeros(3*length(nods),1) ;
% DOFm = zeros(3*length(nodm),1) ;
% for i=1:3
%     % dR = [dR; uB(i,:)'];
%     ROWS = i:3:3*length(nods) ;
%     DOFr(ROWS) = [  3*(nods-1)+i];
%     dR(ROWS) = uB(i,:)' ;
%     ROWS = i:3:3*length(nodm) ;
%     DOFm(ROWS) = [ 3*(nodm-1)+i];
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %
% % 