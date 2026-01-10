function dLref = DistanceMaximPointsEIFEM(Xref,METHOD_COMPUTE_maximum_length)



if METHOD_COMPUTE_maximum_length == 0
    % Method implemented before 13:16
    % -------------------------------
    XrefAUG = [Xref,Xref(:,1)] ;
    dL  =diff(XrefAUG,1,2) ;
    dLref = max(sqrt(sum(dL.^2,1))) ;
else
    % this is actually not needed...
    dLref = 0 ;
    for ipoints  = 1:size(Xref,2)
        distALL = bsxfun(@plus,Xref,Xref(:,ipoints)) ;
        distALLn = max(sqrt(sum(distALL.^2,1))) ;
        dLref = max(dLref,distALLn) ;
    end
    
end
% Compute all the distances
