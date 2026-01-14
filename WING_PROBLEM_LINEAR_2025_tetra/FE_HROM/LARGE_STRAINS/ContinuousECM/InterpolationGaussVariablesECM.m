    function     VARinterp  =  InterpolationGaussVariablesECM(VAR,ECMdata,ngaus,ncomp) 
    % Interpolation of a given variable VAR using the interpolation 
    % information computed with the Continuous Empirical  Cubature Method
    % See 
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
    % ContinuousEmpiricalCubatureM/01_2D_beam_LARGE/README_OECM.pdf
    % JAHO, 15-JAN-2021

    if nargin == 0
        load('tmp1.mat')
    end

    % Set of Gauss points associated to the 
    if ~isempty(ECMdata.setElements)
        setPointsElement = small2large(ECMdata.setElements,ngaus);
        setIndicesLOC = small2large(setPointsElement,ncomp) ;
        VAR = VAR(setIndicesLOC,:);
        %BstZ =  BstRED(setIndicesLOC,:) ;
        nrows =    ncomp*length(ECMdata.setElements) ;
        VARinterp = zeros(nrows,size(VAR,2)) ;
        
        for icomp  =1:ncomp
            INDold = icomp:ncomp:size(VAR,1) ;
            INDnew = icomp:ncomp:nrows ;
            VARinterp(INDnew,:) = ECMdata.Ninterpolation*VAR(INDold,:) ;
        end
        
    else
        VARinterp = [] ; % No need for interpolation 
    end



