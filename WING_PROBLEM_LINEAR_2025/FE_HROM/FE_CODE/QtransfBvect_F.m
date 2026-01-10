function Be = QtransfBvect_F(BeTILDE,ndim)
%
if nargin == 0
    load('tmp2.mat')
end
%Conversion from B-matrix for scalar fields to B-matrix for vector fields
% B matrix such that Fv  =identity +B*d, where Fv stands for the
% deformation gradient
%-------------
% COMPUTATION OF DEFORMATION GRADIENT VECTOR
% ----------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/
% NEW_IDEAS/DOCS/RAUL_BRAVO_THESHIS/ROTATION_SVD_20_Nov_2020/Rotation_SVD.pdf,
% page 33 onwards. See also 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/README_RigidBodyMotions.pdf

nnodeE = size(BeTILDE,2);
nelem = size(BeTILDE,1)/ndim ;
if ndim==2
    nstrain = 4 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
    column1 = 1:2:(nnodeE*2-1) ;
    column2 = 2:2:nnodeE*2 ;
    ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
    ind4 = 4:nstrain:nstrain*nelem ;
    % Fv1 = F11, Fv2 = F22, Fv3 =F12, Fv4 = F21.
    jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; % jnd3 = 3:ndim:ndim*nelem ;
    Be(ind1,column1) = BeTILDE(jnd1,:) ;   % F11
    Be(ind2,column2) = BeTILDE(jnd2,:) ;   % F22
    %  Be(ind3,column1) = BeTILDE(jnd2,:) ;
    Be(ind3,column1) = BeTILDE(jnd2,:) ;   % F12
    Be(ind4,column2) = BeTILDE(jnd1,:) ; % F22
elseif ndim ==3
    
  %  disp('Option not implemented yet')
    
         nstrain = 9 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
         column = cell(1,3) ; 
         column{1} = 1:3:(nnodeE*3-2) ;
         column{2} = 2:3:(nnodeE*3-1) ;
         column{3} = 3:3:nnodeE*3 ;         
         ind = cell(1,9) ;
         for ilocal = 1:length(ind)
             ind{ilocal} = ilocal:nstrain:nstrain*nelem ; 
         end         
         jnd = cell(1,3); 
         for ilocal = 1:length(jnd)
             jnd{ilocal} = ilocal:ndim:ndim*nelem ; 
         end                 
         
         Be(ind{1},column{1}) = BeTILDE(jnd{1},:) ;
         Be(ind{2},column{2}) = BeTILDE(jnd{2},:) ;
         Be(ind{3},column{3}) = BeTILDE(jnd{3},:) ;
         Be(ind{4},column{2}) = BeTILDE(jnd{3},:) ;
         Be(ind{5},column{1}) = BeTILDE(jnd{3},:) ;
         Be(ind{6},column{1}) = BeTILDE(jnd{2},:) ;
         Be(ind{7},column{3}) = BeTILDE(jnd{2},:) ;
         Be(ind{8},column{3}) = BeTILDE(jnd{1},:) ;
         Be(ind{9},column{2}) = BeTILDE(jnd{1},:) ;
     
         
         
    %     ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
    %     ind4 = 4:nstrain:nstrain*nelem; ind5 = 5:nstrain:nstrain*nelem; ind6 = 6:nstrain:nstrain*nelem ;
    %     jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; jnd3 = 3:ndim:ndim*nelem ;
    %     jnd4 = 4:ndim:ndim*nelem; jnd5 = 5:ndim:ndim*nelem; jnd6 = 6:ndim:ndim*nelem ;
    %     Be(ind1,column1) = BeTILDE(jnd1,:) ;
    %     Be(ind2,column2) = BeTILDE(jnd2,:) ;
    %     Be(ind3,column3) = BeTILDE(jnd3,:) ;
    %     Be(ind4,column2) = BeTILDE(jnd3,:) ;
    %     Be(ind4,column3) = BeTILDE(jnd2,:) ;
    %     Be(ind5,column1) = BeTILDE(jnd3,:) ;
    %     Be(ind5,column3) = BeTILDE(jnd1,:) ;
    %     Be(ind6,column1) = BeTILDE(jnd2,:) ;
    %     Be(ind6,column2) = BeTILDE(jnd1,:) ;
else
    error('Incorrect option')
end