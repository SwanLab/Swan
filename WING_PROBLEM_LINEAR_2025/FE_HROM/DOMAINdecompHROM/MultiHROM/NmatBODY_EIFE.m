function [Nmat,w,x,MATPRO,index_elements]  = NmatBODY_EIFE(PhiDEF,PhiRB,PdownsDEF,PdownsRB,CECM_MASS_MATRIX,NstFE,DATAoffline,DATA,MATPRO) ;
% N-matrix for the  empirical interscale FE
% JAHO, 23-Feb-2023
% See comments in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% -------------------------------
if nargin == 0
    load('tmp.mat')
    %  DATAoffline.UseDECMpoints = 1;
end

% Location/weights integration points
% ----------------------------
if  DATAoffline.UseDECMpoints.MASS_MATRIX == 1
    disp('Using DECM for constructing N-matrices')
  %  disp('Error: It must be updated ! ')
    % No continuous ECM (only discrete ECM)
    x = CECM_MASS_MATRIX.xDECM ;
    w = CECM_MASS_MATRIX.wDECM ;
    index_points = CECM_MASS_MATRIX.DECM_indexes_points ;
    index_dofs = small2large(index_points,DATA.MESH.ndim) ;
    Nmat_rb = NstFE(index_dofs,:)*PhiRB*PdownsRB ;
    Nmat_def = NstFE(index_dofs,:)*PhiDEF*PdownsDEF ;
    Nmat = Nmat_rb+Nmat_def;
    index_elements =CECM_MASS_MATRIX.setElements  ; %large2small(index_points,DATA.MESH.ngaus_RHS) ;
    index_strain = small2large(index_points,DATA.MESH.nstrain) ;
    
    
    fff = fieldnames(MATPRO);  % Material properties at the selected integration points
    ngausT_RHS = DATA.MESH.ngausT/DATA.MESH.ngaus_STRESS*DATA.MESH.ngaus_RHS  ;  % 30-May-2025

    for i = 1:length(fff)
        NameProp = fff{i}  ;
        LENGTH_prop = size(MATPRO.(NameProp),1) ;
        if  LENGTH_prop ==  (ngausT_RHS)*DATA.MESH.nstrain
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_strain,:) ;  % DECM
        elseif LENGTH_prop == (ngausT_RHS)
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_points,:) ;
        else
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_elements,:) ;
        end
    end
    
    
    
    
    
else
    
    disp('Using CECM for constructing N-matrices')
    % No continuous ECM (only discrete ECM)
    x = CECM_MASS_MATRIX.xCECM ;
    w = CECM_MASS_MATRIX.wCECM ;
    index_elements = CECM_MASS_MATRIX.setElements ;
    
    
    % This is just for verification purposes
    Nmat_def = NstFE*PhiDEF*PdownsDEF ;
    Nmat_rb = NstFE*PhiRB*PdownsRB ;
    
    Nmat  =  InterpolationGaussVariablesECM(Nmat_def+Nmat_rb,CECM_MASS_MATRIX,DATA.MESH.ngaus_RHS,DATA.MESH.ndim) ;
    %     % Implementation demands a separate representation, see
    %     %/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
    %     Nmat_def_red_all = NstFE*PhiDEF ;
    %     Nmat_def_red  =  InterpolationGaussVariablesECM(Nmat_def_red_all,CECM_MASS_MATRIX,DATA.MESH.ngaus_RHS,DATA.MESH.ndim) ;
    %     Nmat_rb_red_all = NstFE*PhiRB ;
    %     Nmat_rb_red  =  InterpolationGaussVariablesECM(Nmat_rb_red_all,CECM_MASS_MATRIX,DATA.MESH.ngaus_RHS,DATA.MESH.ndim) ;
    
    
    index_points = index_elements*DATA.MESH.ngaus_RHS ;
    index_strain = small2large(index_points,DATA.MESH.nstrain) ;
    
    fff = fieldnames(MATPRO);  % Material properties at the selected integration points
    
    for i = 1:length(fff)
        NameProp = fff{i}  ;
        LENGTH_prop = size(MATPRO.(NameProp),1) ;
        if  LENGTH_prop ==  (DATA.MESH.ngausT)*DATA.MESH.nstrain
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_strain,:) ;  % DECM
        elseif LENGTH_prop == (DATA.MESH.ngausT)
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_points,:) ;
        else
            MATPRO.(NameProp) = MATPRO.(NameProp)(index_elements,:) ;
        end
    end
    
    
    
    
end