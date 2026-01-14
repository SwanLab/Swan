function  DATA_REFMESH = CornersSidesProperties_PLATES(DATA_REFMESH)

if nargin == 0
    load('tmp0.mat')
end
    %% Corner modes (the same for all 4 cornes)
    Vc = ModesCorners_Plates(DATA_REFMESH) ;
    ncorners = 4 ;
    DATA_REFMESH.ModesCornes = cell(ncorners,1) ;
    DATA_REFMESH.ModesCornes(:) = {Vc} ;


    % Side modes  (with respect centroid of each face)
    nface = 4 ;
    DATA_REFMESH.ModesSides = cell(nface,1) ;
    % ----------
    % Sides 1 
    % -------------
    iSIDE = 1;
    [V,As1_1,As1_2] = ModesSides_Plates(iSIDE,DATA_REFMESH) ;
    DATA_REFMESH.ModesSides{1} = V ;
    DATA_REFMESH.ModesSides{3} = V ;
    
     DATA_REFMESH.Acondens_s1 = [As1_1,As1_2] ; 
    
    
    % ----------
    % Sides 2 and 4
    % -------------
    iSIDE = 2;
    [V,As2_1,As2_2 ]= ModesSides_Plates(iSIDE,DATA_REFMESH) ;
    DATA_REFMESH.ModesSides{2} = V ;
    DATA_REFMESH.ModesSides{4} = V ;
    
     DATA_REFMESH.Acondens_s2 = [As2_1,As2_2] ; 

    %% Assembly matrix As 
    ZEROS = zeros(size(As1_1)) ; 

    As = [As1_1   As1_2  ZEROS  ZEROS
          ZEROS   As2_1   As2_2  ZEROS
          ZEROS   ZEROS   As1_2  As1_1
          As2_1   ZEROS   ZEROS  As2_2] ; 

 DATA_REFMESH.Acondens = As ;  
 



%%% BASIS MATRIX THAT RELATES ALL BOUNDARY DOFs with corner amplitudes 
% See pdf document. 
% \begin{equation}
% \BasisINTdom{fc}{e} \defeq \coldos{\BasisINTdom{c}{e}}{\BasisINTdom{s}{e} \Acondens{}{}}


V_c = blkdiag(DATA_REFMESH.ModesCornes{:}) ;  % Modes corners 
V_s = blkdiag(DATA_REFMESH.ModesSides{:}) ; 
% Hence 

V_f = [V_c; V_s*As] ; 

DATA_REFMESH.BasisINTfc = V_f ; 


CHECK_IMPLEMENTATION =0; 

if CHECK_IMPLEMENTATION == 1 
    ndofC = size(DATA_REFMESH.ModesCornes{1},1) ; 
    ndofS1 = size(DATA_REFMESH.ModesSides{1},1) ; 
    ndofS2 = size(DATA_REFMESH.ModesSides{2},1) ; 
    
    nMODES = size(V_f,2) ; 
    
    Vcell = mat2cell(V_f,[ndofC ndofC ndofC ndofC ndofS1 ndofS2 ndofS1 ndofS2],[nMODES]) ; 
    
    
end


