function [MSG,K_axial,K_bending_z,K_bending_y,K_torsion] =...
    BeamPropertiesStiffness(BasisRdef,f1,f2,V,L,AREA,DATA_REFMESH,DATAIN,BasisUdef,K,Tcomp,MSG)

if nargin == 0
    load('tmp1.mat')
end

f = [f1;f2] ;
% -------------------------------------------------------------
KdomRED = BasisUdef'*(K*BasisUdef) ;  % Reduced stiffness matrix
BasisRdef_f = BasisRdef(f,:) ;   % REduced-stiffness at the boundary
BasisUdef_f = BasisUdef(f,:) ;
%  \Hqr{e} =  \BasisUdefT{e}{f}  \BasisRdef{e}{f}
Hqr = BasisUdef_f'*BasisRdef_f ;




%  \Kbeam{e}{} \defeq (\Hqr{e^T}\KdomRED{e^{-1}}{} \Hqr{e})^{-1},
Kbeam_inv = Hqr'*(KdomRED\Hqr) ;
Kbeam = inv(Kbeam_inv) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \BasisINTdom{1 }{e^T} (\BasisRdef{e}{\face{1}}\Kbeam{e}{}\BasisRdef{e^T}{\face{1}}) \BasisINTdom{1}{e}
iface= 1;
VB1 = V{iface}'*BasisRdef(f1,:) ;
iface= 2;
VB2 = V{iface}'*BasisRdef(f2,:) ;
Kskel_11 = VB1*Kbeam*VB1' ;
Kskel_12 = VB1*Kbeam*VB2' ; ;
Kskel_21 = VB2*Kbeam*VB1' ; ;
Kskel_22 =VB2*Kbeam*VB2' ; ;


%%%% STIFNESS MATRIX % ------------------------------
Kskel = [Kskel_11, Kskel_12; Kskel_21, Kskel_22] ;

MSG{end+1} = '----------------------------------';
MSG{end+1} = 'GEOMETRIC PROPERTIES SLICE ';
MSG{end+1} = '----------------------------------';
MSG{end+1} =['LENGTH =',num2str(L)] ;

%DATAIN = DefaultField(DATAIN,'MASS_MATRIX_FORMULATION_REACTIONS',0) ;
%%Removed 18-th-Dec-2019
ndim = size(DATA_REFMESH.COOR,2) ;

if  ~isempty(DATA_REFMESH.GEOMETRIC_PROPERTIES_VOLUME)   % && DATAIN.MASS_MATRIX_FORMULATION_REACTIONS == 1
    I = DATA_REFMESH.GEOMETRIC_PROPERTIES_VOLUME  ;
    if ndim == 3
        Inert = I(4:6,4:6) ;
        MSG{end+1} =['VOLUME =',num2str(I(1,1))] ;
        MSG{end+1} ='MATRIX OF INERTIAS ';
        MSG{end+1} ='------------------------------' ;
        MSG{end+1} =num2str(Inert,3)  ;
    else
        MSG{end+1} =['Area (slice) = ',num2str(I(1,1))]; ;
        MSG{end+1} =['Moment of Inertia (slice) = ',num2str(I(3,3))]  ;
        
    end
    
end
MSG{end+1}='----------------------------------';

MSG{end+1}='----------------------------------';
MSG{end+1}='GEOMETRIC PROPERTIES CROSS-SECTION  ';
MSG{end+1}='----------------------------------';
MSG{end+1}='-----------------------------';
MSG{end+1}=['AREA1 =',num2str(DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{1}(1))];
MSG{end+1}=['AREA2 =',num2str(DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{2}(1))];
MSG{end+1}='-----------------------------';
if size(DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{1},1) == 6
    COLS = [4:6] ;
else
    COLS = [3] ;
end
INERT1 = DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{1}(COLS,COLS) ;
INERT2= DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{2}(COLS,COLS) ;

MSG{end+1}='MATRIX OF INERTIAS, face A ';
MSG{end+1}='------------------------------';
MSG{end+1}=num2str(INERT1,4) ;
MSG{end+1}=['MATRIX OF INERTIAS, face B'];
MSG{end+1}='------------------------------';
MSG{end+1}=num2str(INERT2,4);

SHOW_K = 1;
if SHOW_K ==1
    MSG{end+1}='BEAM STIFFNESS MATRIX';
    MSG{end+1}='K_11';
    MSG{end+1}=num2str(Kskel_11,6);
    MSG{end+1}='K_12' ;
    MSG{end+1}=num2str(Kskel_12,6);
    MSG{end+1}='K_21' ;
    MSG{end+1}=num2str(Kskel_21,6);
    MSG{end+1}='K_22' ;
    MSG{end+1}=num2str(Kskel_22,6);
    
end


if abs(DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{1}(1) - DATA_REFMESH.GEOMETRIC_PROPERTIES_CROSS_SECTION{2}(1)) <1e-3
    
    
    if size(Kskel_11,1) == 6
        [K_axial,K_bending_z,K_bending_y,K_torsion,MSG] = PropStiffBeam3D(Kskel_11,Kskel_12,Kskel,L,MSG) ;
        
    elseif size(Kskel_11,1) == 3
        [K_axial,K_bending_z,MSG ] = PropStiffBeam2D(Kskel_11,Kskel_12,Kskel,L,MSG) ;
        
    end
    
end
