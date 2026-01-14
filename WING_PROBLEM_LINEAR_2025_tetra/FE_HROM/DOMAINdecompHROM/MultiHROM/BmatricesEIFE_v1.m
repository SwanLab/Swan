function [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED]  = ...
    BmatricesEIFE_v1(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH) 
% B-matrix empirical interscale FE
% JAHO, 23-Feb-2023/1-MArch-2023
% Updated, improved version of BmatricesEIFE.m 
% See comments in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% -------------------------------
if nargin == 0
    load('tmp1.mat')
 %   DATAoffline.UseDECMpoints = 0;
end

% Location/weights integration points
% ----------------------------
%DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',0) ;

if ~isempty(CECM_intforces.xCECM)  && DATAoffline.UseDECMpoints == 1
    disp('Using DECM for constructing B-matrices')
    % No continuous ECM (only discrete ECM)
    x = CECM_intforces.xDECM ;
    w = CECM_intforces.wDECM ;
    index_points = CECM_intforces.DECM_indexes_points ;
    index_strain = small2large(index_points,DATA.MESH.nstrain) ;
    Bmat = BstFE(index_strain,:)*PhiDEF*PdownsDEF ;
     BmatRED = BstFE(index_strain,:)*PhiDEF; % *PdownsDEF ;
    index_elements = large2small(index_points,DATA.MESH.ngaus_STRESS) ;
    
    
     % RECONSTRUCTION OF PK2 STRESSES
    % Stress basis matrix at the selected points
    BasisStwo_z = BasisStwo(index_strain,:)  ; 
    
else
    
    
    disp('Using CECM for constructing B-matrices')
    % No continuous ECM (only discrete ECM)
    x = CECM_intforces.xCECM ;
    w = CECM_intforces.wCECM ;
    index_elements = CECM_intforces.setElements ;
    index_points = index_elements*DATA.MESH.ngaus_STRESS ;
    index_strain = small2large(index_points,DATA.MESH.nstrain) ;
    
    %% THIS IS THE B-matrix of the parent element. It cannot be used directly for elements not oriented like
    % the parent element. Nevertheless, we compute for verification
    % purposes
    Bmat_all = BstFE*PhiDEF*PdownsDEF ;
    Bmat  =  InterpolationGaussVariablesECM(Bmat_all,CECM_intforces,DATA.MESH.ngaus_STRESS,DATA.MESH.nstrain) ;
    
    % GENERAL IMPLEMENTATION
    % We begin by defining the reduced Bmatrix
    BmatRED_all =  BstFE*PhiDEF ;  
    % Now we interpolate this matrix 
    BmatRED  =  InterpolationGaussVariablesECM(BmatRED_all,CECM_intforces,DATA.MESH.ngaus_STRESS,DATA.MESH.nstrain) ;
% 
%     
%     % Separated representation (pending vectorization) 
%     BePOINT_all = [] ; 
%     PhiDEFelem = cell(length(w),1) ; 
%     % ------------------------
%     for ipoint = 1:length(w)
%         ielem = index_elements(ipoint) ;
%         nodesELEM = MESH.CN(ielem,:); 
%         DOFsELEM  = small2large(nodesELEM,DATA.MESH.ndim) ; 
%         gaussELEM = small2large(ielem,DATA.MESH.ngaus_STRESS) ; 
%         istrainELEM = small2large(gaussELEM,DATA.MESH.nstrain) ; 
%         BeELEM = BstFE(istrainELEM,DOFsELEM); 
%         BePOINT_all = [BePOINT_all ; BeELEM] ; 
%         PhiDEFelem{ipoint} = PhiDEF(DOFsELEM,:); 
%     end
%     nrows = DATA.MESH.nstrain*length(w) ; 
%     BmatFE = zeros(nrows,size(BePOINT_all,2)) ;
%     ncomp = DATA.MESH.nstrain ; 
%     for icomp  =1:ncomp
%         INDold = icomp:ncomp:size(BePOINT_all,1) ;
%         INDnew = icomp:ncomp:nrows ;
%         BmatFE(INDnew,:) = CECM_intforces.Ninterpolation*BePOINT_all(INDold,:) ;
%     end
%     rows = DATA.MESH.nstrain*ones(1,length(w)) ;
%     cols = [1,size(BmatFE,2)] ; 
%     BmatFE = mat2cell(BmatFE,rows',size(BmatFE,2)) ; 
%     Bgrad = cell(size(BmatFE)) ; 
    
%    % BmatFE_Phi = cell(size(Bgrad)) ; 
%     for  ipoint = 1:length(w)
%         
%     %    BmatFE_Phi{ipoint} = BmatFE{ipoint}*PhiDEFelem{ipoint} ; 
%         
%         Bgrad{ipoint} = zeros(DATA.MESH.ndim,size(BmatFE{ipoint},2)/DATA.MESH.ndim) ; 
%          
%         for idim = 1:DATA.MESH.ndim
%             Bgrad{ipoint}(idim,:) = BmatFE{ipoint}(idim,idim:DATA.MESH.ndim:end) ; 
%         end
%           
%     end 
%     
    
    
    
    
    
    % RECONSTRUCTION OF PK2 STRESSES
    % Stress basis matrix at the selected points
    BasisStwo_z  =  InterpolationGaussVariablesECM(BasisStwo,CECM_intforces,DATA.MESH.ngaus_STRESS,DATA.MESH.nstrain) ;
    
    
    %     fff = fieldnames(MATPRO);  % Material properties at the selected integration points
    %
    %     for i = 1:length(fff)
    %         NameProp = fff{i}  ;
    %         LENGTH_prop = size(MATPRO.(NameProp),1) ;
    %         if  LENGTH_prop ==  (DATA.MESH.ngausT)*DATA.MESH.nstrain
    %             MATPRO.(NameProp) = MATPRO.(NameProp)(index_strain,:) ;  % DECM
    %         elseif LENGTH_prop == (DATA.MESH.ngausT)
    %             MATPRO.(NameProp) = MATPRO.(NameProp)(index_points,:) ;
    %         else
    %             MATPRO.(NameProp) = MATPRO.(NameProp)(index_elements,:) ;
    %         end
    %     end
    %
    
    
end


RECONS_PK2stresses = [];
RECONS_PK2stresses.coeff  =  (BasisStwo_z'*BasisStwo_z)\BasisStwo_z' ;
RECONS_PK2stresses.BASIS = BasisStwo ;

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