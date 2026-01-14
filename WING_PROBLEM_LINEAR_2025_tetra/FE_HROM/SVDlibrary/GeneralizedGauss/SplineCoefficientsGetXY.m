function [S,Ux_splineCOEFF,Ux_splineCOEFF_der,Vy_splineCOEFF,Vy_splineCOEFF_der,Ux,Vy] =...
    SplineCoefficientsGetXY(x,y,Gi,TOLfitting)


ly = length(y) ;
    C = reshape(Gi,ly,[])';      
    % ---------------------------------------
    % SVD 
    
    DATASVD.RELATIVE_SVD = 1;
    [Ux,S,Vy] = RSVDT(C,TOLfitting,[],0,DATASVD) ;
    
%     % Avoid changes of order of magnitude in the SVD 
%     DATALOC = DefaultField(DATALOC,'TruncationCriterionSVDmodesJumpOrderMagnitudeLOG_fitting',5) ;
%     if ~isempty(DATALOC.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG_fitting)
%         Sincre = log10(S(1:end-1)./S(2:end)) ;
%         III = find(Sincre >= DATALOC.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG_fitting) ;
%         if ~isempty(III)
%             nmodes= III(1) ;
%             Ux = Ux(:,1:nmodes) ;
%             S = S(1:nmodes) ;
%             Vy = Vy(:,1:nmodes) ;
%         end
%     end  
    % ----------------------------------------------------------------------------------------
    
    % Spline coefficients 
    
    Ux_splineCOEFF = cell(length(S),1) ;
    Vy_splineCOEFF = cell(length(S),1) ;  
    
    Ux_splineCOEFF_der = cell(length(S),1) ;
    Vy_splineCOEFF_der = cell(length(S),1) ;
     
    for i = 1:length(S)
        Ux_splineCOEFF{i} =  spline(x,Ux(:,i)) ;
        Vy_splineCOEFF{i} = spline(y,Vy(:,i)) ;
        Ux_splineCOEFF_der{i} = ppDer(Ux_splineCOEFF{i}) ;
        Vy_splineCOEFF_der{i} = ppDer(  Vy_splineCOEFF{i});
        
    end