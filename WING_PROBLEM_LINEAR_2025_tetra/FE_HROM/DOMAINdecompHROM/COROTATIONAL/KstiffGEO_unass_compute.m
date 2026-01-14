function [KcGEOunassemGLOloc,PdownsRBcoupROTc_i ]= KstiffGEO_unass_compute(OPERFE,DATA,D_QrotALL,FintCunassembLOC,Qrot)
%  Unassembled geometric stiffness matrix (GLO-LOC)
% See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 19-feb-2025, WEDNESDAY, 11:16, Balmes 185,  Barcelona/25-fEB-2025, TUESDAY,UPC TERRASSA.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/06_COROT_3D.mlx
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
if nargin == 0
    load('tmp1.mat')
end


FintCunassemb  =  D_QrotALL*FintCunassembLOC ;
D_FintCunassembGLO = DiagC_FintCunassemb(FintCunassemb,OPERFE.INDEXsparseFINT)  ;
if DATA.MESH.ndim ==2
    %      \KcGEOunassemGLOloc \defeq   \DiagC{\AspMAT}  \DiagC{\FintCunassemb}  \DiagC{\PdownsRBlROT}
    KcGEOunassemGLOloc =  OPERFE.D_AspMAT*(D_FintCunassembGLO*OPERFE.D_PdownsRBlROT)  ;
    PdownsRBcoupROTc_i = [] ;
else
    % ------------
    % 3D PROBLEMS
    %-------------
    %\{\KcGEOunassemGLOloc =  \sum_{i=1}^3  \Par{ \DiagC{\AspMATc{i}}   \DiagC{\FintCunassemb}    \DiagC{\PdownsRBcoupROTc{i}}  }
    %   \DiagC{\PdownsRBcoupROTc{i}} = \DiagC{[\Qrot]_i}   \DiagC{\PdownsRBlROT},  \hspace{0.5cm}  i  =1,2,3
    
    % Computing \DiagC{[\Qrot]_i}
    D_Qrot_i = Qrot_by_rows_element(Qrot,OPERFE.INDEXsparseROTelem_byROWS) ;
    ndim = 3;
    PdownsRBcoupROTc_i = cell(size(D_Qrot_i)) ;
    for idim = 1:ndim
        PdownsRBcoupROTc_i{idim} =  D_Qrot_i{idim}*OPERFE.D_PdownsRBlROT ;
        if idim == 1
            KcGEOunassemGLOloc = OPERFE.D_AspMAT{idim}*(D_FintCunassembGLO*PdownsRBcoupROTc_i{idim}) ;
        else
            KcGEOunassemGLOloc = KcGEOunassemGLOloc +  OPERFE.D_AspMAT{idim}*(D_FintCunassembGLO*PdownsRBcoupROTc_i{idim}) ;
        end
    end
    
end
