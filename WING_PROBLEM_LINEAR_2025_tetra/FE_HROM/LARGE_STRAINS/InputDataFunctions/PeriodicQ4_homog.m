function  [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] =...
    PeriodicQ4_homog(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,DOFr_LINEAR,dR_LINEAR)
% Adaptation  PeriodicQ4 for homogenization problems
% JAHO, 7-Oct-2025, UPC, Terrassa.  

%
%  PeriodicQ4   returngs the matrices and vectors defining affine boundary
% conditions for EIF Q4 elements with corners
% JAHO, 17-Nov-2023, CIMNE, UPC (Barcelona)
% See details in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
% ------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    
end

 COMB_FACES= {[1,2],[2,3],[3,4],[4,1]} ;  % Faces defining corners (if any)
 % CHECK IF THERE ARE CORNERS
icorner = 1;
ind_FACE_A = COMB_FACES{icorner}(1) ;
ind_FACE_B = COMB_FACES{icorner}(2) ;
nodes_FACE_A = MESH.NODES_FACES{ind_FACE_A}  ;
nodes_FACE_B = MESH.NODES_FACES{ind_FACE_B}  ;
III = intersect(nodes_FACE_A,nodes_FACE_B) ;

%DATALOC = DefaultField(DATALOC,'Ufluct',[])  ; 

% if isempty(DATALOC.Ufluct)
%     
%     if ~isempty(III) || DATALOC.PRESCRIBED_DIRECTLY_FLUCTUATIONS_AT_END_POINTS_WHEN_NO_CORNERS == 1
%         % CASE WITH CORNERS (NODES OF THE A4 ELEMENTS ARE MATERIAL POINTS )

DATALOC.PERIODICITY_FACES = {[1,3],[2,4]} ; 
        [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4_cornersHOM(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,COMB_FACES,...
            DOFr_LINEAR,dR_LINEAR)  ;
%     else
%         disp(['It has been proven that this option is equivalent to PeriodicQ4_corners...'])
%         [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4_NOcorners(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,COMB_FACES);
%         
%     end
%     
% else
% % Prescribed fluctuations 
%     [DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PrescFluctQ4(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,COMB_FACES)  ;
%     
% end
%  