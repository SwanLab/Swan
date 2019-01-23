function [d_u,free1] = generalized_solver(d_u,StifMat,element,free_df,fix_df,dim,fixnodes,fextlod)

if (size(element.pnods,2)>0 || size(element.lglib,2)>0)
    % periodic conditions
    [d_u,free1] = solverp(d_u,StifMat,element.pnods,fextlod,fix_df,dim);
else
    d_u = solver(d_u,free_df,StifMat,fextlod,fix_df,fixnodes);
    free1 = [];
end
    
end