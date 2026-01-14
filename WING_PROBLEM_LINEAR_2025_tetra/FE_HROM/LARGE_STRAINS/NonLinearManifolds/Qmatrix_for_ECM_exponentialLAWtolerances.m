function [Q,S,V] =  Qmatrix_for_ECM_exponentialLAWtolerances(SNAPredFINT_nw,sqrt_wST,DATAoffline)

    SNAPredFINT = cell(size(SNAPredFINT_nw));
    rNORM_snapFint = zeros(size(SNAPredFINT_nw)) ;
    for iproj = 1:length(SNAPredFINT_nw)
        SNAPredFINT{iproj} = bsxfun(@times,SNAPredFINT_nw{iproj},sqrt_wST) ;
        rNORM_snapFint(iproj) = norm(SNAPredFINT{iproj},'fro')^2 ;
    end
    
    rMAX= max(rNORM_snapFint)  ;
    rNORM_snapFint = sqrt(rNORM_snapFint/rMAX) ;
    % The above is the ratio between the norm of a particular block and the
    % global norm
    % This variable will be used in adjusting local tolerances, as explained in
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/11_MAW_ECM_plast.mlx
    % and
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/Tolearnces_SRSVD.pdf
    
    % Proposed law
    p = DATAoffline.Exponent_Function_relating_global_local_TOL_fint ;
    TOL = DATAoffline.errorFINT ;
    TOL_i = 1 + (TOL-1)*rNORM_snapFint.^p ;
    PLOT_TOL = 1;
    if  PLOT_TOL == 1
        figure(3467)
        hold on
        xlabel('Snapshot')
        ylabel('Local Tolerance')
        plot(TOL_i)
        title(['Local tolerances used for each block (internal force matrix)'])
    end
    
    %TOL_loc = ones(length(SNAPredFINT),1)*DATAoffline.errorFINT ;
    DATAsvdLOC.HIDE_OUTPUT =  1;
    [Q,S,V] = SRSVD(SNAPredFINT,TOL_i,DATAsvdLOC) ;