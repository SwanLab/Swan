function [SNAPbasic,SNAPcompl,DATAoffline] = ScalingSnapshotsDISP(SNAPdisp,SNAPdisp_AUX,INDEX_BASIC_TESTS,DATAoffline,INDEX_COMPL_TESTS)
%--------------------------------------------------------------------------
% FUNCTION: ScalingSnapshotsDISP
%
% PURPOSE:
%   Scales displacement snapshot data (e.g., from training simulations) 
%   according to predefined influence weights for each domain, allowing 
%   the user to modulate the contribution of different domains in the 
%   computation of deformation modes (e.g., in EIFEM or ROM settings).
%
%   The main domain (often considered the reference) and any auxiliary 
%   domains are scaled independently using user-defined or default 
%   scaling coefficients (between 0 and 1, default = 1). This supports
%   weighted multi-domain training and enhanced control over basis
%   construction.
%
% INPUTS:
%   - SNAPdisp           : Cell array of displacement snapshots, 
%                          combining all (main + auxiliary) domains.
%   - INDEX_BASIC_TESTS  : Vector of indices corresponding to elastic tests.
%   - DATAoffline        : Structure containing:
%         • .INFLUENCE_SNAPSHOTS_DISP.FIRST_DOMAIN : Scaling factor for main domain (default = 1)
%         • .INFLUENCE_SNAPSHOTS_DISP.AUXDOMAINS   : Vector of scaling factors for each auxiliary domain
%
% GLOBAL VARIABLES (must be available in caller workspace or defined prior):
%   - INDEX_COMPL_TESTS : Indices for inelastic tests in SNAPdisp
%   - SNAPdisp_AUX      : Cell array of snapshots from auxiliary domains (optional)
%
% OUTPUTS:
%   - SNAPbasic : Scaled and concatenated basic (elastic) displacement snapshots
%   - SNAPcompl : Scaled and concatenated complementary (inelastic) snapshots
%
% REMARKS:
%   - Elastic snapshots from auxiliary domains are stored as 
%     `SNAPdisp_AUX(:,INDEX_BASIC_TESTS)`, and similarly for complementary.
%   - If no auxiliary data is present, only the main domain’s snapshots 
%     are returned.
%   - Scaling factors are useful when combining simulations with different
%     boundary conditions, resolutions, or material models.
%
% EXAMPLE USAGE:
%   [SNAPbasic, SNAPcompl] = ScalingSnapshotsDISP(SNAPdisp, INDEX_BASIC_TESTS, DATAoffline);
%
% RELATED FILES:
%   - /109_EIFEM_largeROT/13_auxALEX.mlx
%   - /109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   First version: 1-June-2025
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end



SNAPbasic = cell2mat(SNAPdisp(INDEX_BASIC_TESTS)) ;  % ELASTIC MODES

 
%--------------------------------------------------------------------------
DATAoffline = DefaultField(DATAoffline,'INFLUENCE_SNAPSHOTS_DISP',[]);
DATAoffline.INFLUENCE_SNAPSHOTS_DISP = DefaultField(DATAoffline.INFLUENCE_SNAPSHOTS_DISP,'FIRST_DOMAIN',1);
SNAPbasic = SNAPbasic*DATAoffline.INFLUENCE_SNAPSHOTS_DISP.FIRST_DOMAIN ;

SNAPcompl = SNAPdisp(INDEX_COMPL_TESTS) ;
for isnap = 1:length(SNAPcompl)
    SNAPcompl{isnap} = SNAPcompl{isnap}*DATAoffline.INFLUENCE_SNAPSHOTS_DISP.FIRST_DOMAIN ;    
end


if ~isempty(SNAPdisp_AUX)
    % ADDITIONAL DOMAINS COMES INTO PLAY
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
    DATAoffline.INFLUENCE_SNAPSHOTS_DISP = DefaultField(DATAoffline.INFLUENCE_SNAPSHOTS_DISP,'AUXDOMAINS',ones(size(SNAPdisp_AUX)));
    % Influence of each domain in the the final set of modes (scaling factor between 0 and 1, default 1)
    SNAPbasic_aux =  SNAPdisp_AUX(:,INDEX_BASIC_TESTS)  ;    
    for iprojAUX = 1:size(SNAPbasic_aux,2)
        for iaux=  1:size(SNAPbasic_aux,1)
            SNAPbasic_aux{iaux,iprojAUX} = SNAPbasic_aux{iaux,iprojAUX}*DATAoffline.INFLUENCE_SNAPSHOTS_DISP.AUXDOMAINS(iaux) ;
        end
    end    
    SNAPbasic = [SNAPbasic,cell2mat(SNAPbasic_aux(:)')] ;
    if ~isempty(SNAPcompl)
        SNAPcompl_aux =  SNAPdisp_AUX(:,INDEX_COMPL_TESTS)  ;        
        for iprojAUX = 1:size(SNAPcompl_aux,2)
            for iaux=  1:size(SNAPcompl_aux,1)
                SNAPcompl_aux{iaux,iprojAUX} = SNAPcompl_aux{iaux,iprojAUX}*DATAoffline.INFLUENCE_SNAPSHOTS_DISP.AUXDOMAINS(iaux) ;
            end
            SNAPcompl{iprojAUX} = [SNAPcompl{iprojAUX},cell2mat(SNAPcompl_aux(:,iprojAUX)')] ; 
        end               
      %  SNAPcompl = [SNAPcompl,cell2mat(SNAPcompl_aux(:)')] ;
    end
    
end
