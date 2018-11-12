function projection_Ch = projectionCh(inv_matCh,alpha,beta)

weights = alpha*beta';
projection_Ch = sum(sum(weights.*inv_matCh));
    
end