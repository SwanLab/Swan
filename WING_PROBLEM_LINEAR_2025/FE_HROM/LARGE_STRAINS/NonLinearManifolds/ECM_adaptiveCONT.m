function ECMdata = ECM_adaptiveCONT(DATAoffline,DATA_GENGAUSS,A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE)
%--------------------------------------------------------------------------
% function ECMdata = ECM_adaptiveCONT(DATAoffline,DATA_GENGAUSS,...
%     A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE)
%
% PURPOSE:
%   Driver routine for constructing an adaptive Empirical Cubature Method
%   (ECM) rule tailored to the given internal force snapshots. The function
%   selects a reduced set of quadrature points and weights for hyperreduced
%   finite element simulations. Depending on the offline options, either a
%   standard point-based ECM is applied or (legacy) element-based and
%   continuous ECM branches are invoked.
%
% INPUTS:
%   DATAoffline   : structure with offline options
%       • USE_ELEMENT_BASED_APPROACH – switch between point- and
%         element-based ECM (default: 0 = point-based).
%       • Hyperreduction_METHOD_projection – string flag to control
%         manifold-aware ECM projection ("MANIFOLD") or standard ECM.
%   DATA_GENGAUSS : structure controlling generalized Gauss integration
%       • ACTIVE – if nonzero, continuous ECM branch is attempted
%         (currently disabled).
%   A_internalFORCES_ECM_lin : matrix of linear internal force snapshots.
%   A_internalFORCES_ECM_non : matrix of nonlinear internal force snapshots.
%   DATA          : global model/mesh information (including mesh and gauss
%                   point counts).
%   OPERFE        : structure with FE operators (e.g. quadrature weights).
%
% OUTPUTS:
%   ECMdata : structure with all information about the ECM rule
%       • setPoints       – indices of selected quadrature points.
%       • wRED            – reduced positive weights.
%       • setElements     – parent element IDs of selected points.
%       • setPoints_slv   – (if MANIFOLD) slave quadrature points.
%       • wRED_slv        – (if MANIFOLD) slave weights.
%       • setElements_slv – (if MANIFOLD) parent elements of slave points.
%       • DATA_regress_eta_der – regression structure for nonlinear
%                                master/slave mapping (if MANIFOLD).
%       • etaNON          – empty if no manifold mapping is used.
%
% NOTES:
%   • The default and maintained branch is the point-based Discrete ECM.
%   • The "MANIFOLD" projection option attempts a nonlinear mapping
%     between master and slave ECM points. This is deprecated as of
%     2-SEP-2025 due to poor performance in HROM applications.
%   • Continuous and element-based ECM branches are preserved as error
%     stubs for reference but are no longer maintained (since 2-JUL-2025).
%
% HISTORY:
%   Written as adaptive extension of Discrete ECM with optional manifold
%   master/slave mapping. 
%   JAHO, Cartagena–Barcelona, July–September 2025
%--------------------------------------------------------------------------

DATAoffline = DefaultField(DATAoffline,'USE_ELEMENT_BASED_APPROACH',0);

% For standard usage we do point‑based Discrete ECM
if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
    if DATA_GENGAUSS.ACTIVE == 0
        % ----- Discrete ECM over A_fint -----
        [setPoints,wRED,ERROR_GLO,DATAOUT] = ...
            DiscreteECM_givenAmat_nonFINT(A_internalFORCES_ECM_lin,A_internalFORCES_ECM_non,DATA,OPERFE.wSTs,DATAoffline);
        
        ECMdata.setPoints = setPoints;
        ECMdata.wRED      = wRED;
        
        proporP = length(setPoints)/DATA.MESH.ngausT*100;
        disp(['Number of ECM points = ',num2str(length(setPoints)), ...
            ' (',num2str(proporP),' % total)'])
        
        % Also report selected elements for convenience
        setElements = large2smallREP(setPoints,DATA.MESH.ngaus);
        disp('****************************+')
        disp(['List of selected m = ',num2str(length(setElements)),' elements'])
        disp(num2str(setElements'))
        clipboard('copy',num2str(setElements'));              % quick paste into notes
        ECMdata.setElements = setElements;
        
        % ----- Optional manifold‑aware ECM (master/slave mapping) -----
        switch DATAoffline.Hyperreduction_METHOD_projection
            case 'MANIFOLD'
                warning('Option no longer maintained... It proved not viable at the HROM stage (2-09-2025)')
                disp('Master/Slave nonlinear mapping between ECM points')                
                %
                [ECMdata.setPoints,ECMdata.setPoints_slv, ...
                    ECMdata.wRED,ECMdata.wRED_slv, ...
                  ECMdata.DATA_regress_eta_der] = ...
                    DiscreteECM_adaptWEIGHTS_ELPLAST2(A_internalFORCES_ECM_non, ...
                    setPoints,wRED, ...
                    DATA_interp,OPERFE.wSTs); 
                
                ECMdata.setElements     = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus);
                disp(['Master element(s) = ',num2str(ECMdata.setElements(:)')])
                ECMdata.setElements_slv = large2smallREP(ECMdata.setPoints_slv,DATA.MESH.ngaus);
                disp(['Slave element(s)  = ',num2str(ECMdata.setElements_slv(:)')])
            otherwise
                ECMdata.etaNON = []; % STANDARD projection: no manifold mapping
        end
    else
        % Continuous ECM branch – kept as placeholder (not maintained)
        error('Continuous ECM option not available (2‑JUL‑2025)')
    end
else
    % Element‑based approach (legacy; not maintained)
    error('Element‑based ECM option not maintained (2‑JUL‑2025)')
end