function [Fb Nst  DATA wSTs posgp]= ComputeFbVect(COOR,CN,TypeElement,fNOD,DATA,Nst,wSTs,posgp)
%%%%
% This subroutine   returns the  body force    contribution (Fb)
% due to the % global external force vector. Inputs:   COOR, CN,
%TypeElement,
% fNOD (nnode*ndim x 1):   Body force function at the nodes  of the mesh.
% wST : Vector of weights
%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 27-Oct-2015
%dbstop('11')
if nargin == 0
    load('tmp1.mat')
end
%%% Dimensions of the problem
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;




Nst = [] ;

% 21-Apr-2016, Number of Gauss points for integrating RHS and mass MATRIX may be different
% from number of Gauss points used in integrating stiffness matrix


DATA = DefaultField(DATA,'posgp_given',[]) ; 
% 4-Nov-2022 --- GAuss points given provided by the user (see for instance
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/linRECT_CROSS_SECTION/DATAFE_slice_linear3D.m
%
if isempty(DATA.posgp_given)
TypeIntegrand = 'RHS' ;
else
   TypeIntegrand = {DATA.posgp_given,DATA.weights_given} ;  
end

[  wSTs wST posgp] = ComputeW_RHS(COOR,CN,TypeElement,ndim,TypeIntegrand)  ;

DATA = DefaultField(DATA,'CalculateNst',1) ;

if any(fNOD)   || DATA.CALCULATE_MASSMATRIX  == 1 || DATA.CalculateNst == 1
    %------------------------------------------------------------
    % Matrix consisting of all element shape function matrices
    % -----------------------------------------------------------
    if isempty(Nst)
        disp('Computing N-matrices for all elements...')
        [ Nelem,posgp  ] = ComputeNelemALL(TypeElement,nnodeE,ndim,nelem,TypeIntegrand) ;
        % -------------------------------------------------------------------------
        disp('Done')
    end
    % Assembly of Nst
    %   dbstop('32')
    ngaus = size(posgp,2) ;
    DATA.ngaus = ngaus ;
    if isempty(Nst)
    disp('Assembly of Nst (stacked shape functions matrix)...')
    Nst = AssemblyNGlobal(Nelem,nelem,nnodeE,ndim,ngaus,CN,nnode) ;
    end
    disp('Done')
    % Nst including weights
    % Matrix of weights
    disp('COmputing NstW')
    wDIAG = CompWeightDiag(wSTs,ndim)  ;
    NstW = wDIAG*Nst ;
    disp('Done')
    %--------------------------------------------------------------------------
    % - Computing contribution to Fb of external nodal forces defined on nodes
    % In a structural problem, these forces are normally zero. The only
    % volumetric force is self-weigth, which is calculated after the below
    % operations
    % -------------------------------------------------------------------------
    if  any(fNOD)
        disp( 'Computing body forces')
        NF = Nst*fNOD ;
        Fb = NstW'*NF ;
    else
        Fb = zeros(ndim*nnode,1) ; 
    end
    %%%%%% SELF-WEIGHT %%%%%%%%%%%%%%
    DATA = DefaultField(DATA,'INCLUDE_SELFWEIGTH',0) ;
    if DATA.INCLUDE_SELFWEIGTH == 1
        % Global density vector (defined for all elements)
        densGLO = DATA.densGLO ;
        if isempty(densGLO)
            disp(['Self-weight cannot be computed because there are no density data'])
            disp(['Define the density of each material in input variable: MATERIAL.PLY(imaterial).densCOMP  '])
            error([''])
        end
        % ----------------------------
        densGAUSS = zeros(ngaus*length(densGLO),1) ;
        for igaus = 1:ngaus
            densGAUSS(igaus:ngaus:end) = densGLO ;
        end
        %         Fself = NstW'*
        DATA = DefaultField(DATA,'gravity_acceleration',[0 -9.81 0]) ; % m/s^2
        FORCE_GRAVITY = zeros(ngaus*ndim*nelem,1) ;
        for idim = 1:ndim
            FORCE_GRAVITY(idim:ndim:end) = DATA.gravity_acceleration(idim)*densGAUSS ;
        end
        F_SELFWEIGHT = NstW'*FORCE_GRAVITY ;
        
            Fb = Fb + F_SELFWEIGHT ;

    end
    
    
    
    
    
else
    Fb = zeros(ndim*nnode,1) ;
end








