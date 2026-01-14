function [ECMdata,setCandidates,LIST_OF_CANDIDATES] = SAW_ECM_elastplast(DATAoffline,DATA_GENGAUSS,A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE,qPLAST)
 %---------------------
 if nargin == 0
     load('tmp2.mat')
 end

DATAoffline = DefaultField(DATAoffline,'USE_ELEMENT_BASED_APPROACH',0);

% For standard usage we do point‑based Discrete ECM
if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
    if DATA_GENGAUSS.ACTIVE == 0
        % ----- Discrete ECM over A_fint -----
      % [ECMdata,setCandidates,LIST_OF_CANDIDATES] = ...
      %      SAW_ECM_elastplastLOC(A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE.wSTs,DATAoffline,qPLAST);
        
         ECMdata= ...
            SAW_ECM_elastplastLOC2(A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE.wSTs,DATAoffline,qPLAST);
        
%         ECMdata.setPoints = setPoints;
%         ECMdata.wRED      = wRED;
%         
%         proporP = length(setPoints)/DATA.MESH.ngausT*100;
%         disp(['Number of ECM points = ',num2str(length(setPoints)), ...
%             ' (',num2str(proporP),' % total)'])
%         
%         % Also report selected elements for convenience
%         setElements = large2smallREP(setPoints,DATA.MESH.ngaus);
%         disp('****************************+')
%         disp(['List of selected m = ',num2str(length(setElements)),' elements'])
%         disp(num2str(setElements'))
%         clipboard('copy',num2str(setElements'));              % quick paste into notes
%         ECMdata.setElements = setElements;
        
        
         
    else
        % Continuous ECM branch – kept as placeholder (not maintained)
        error('Continuous ECM option not available (2‑JUL‑2025)')
    end
else
    % Element‑based approach (legacy; not maintained)
    error('Element‑based ECM option not maintained (2‑JUL‑2025)')
end