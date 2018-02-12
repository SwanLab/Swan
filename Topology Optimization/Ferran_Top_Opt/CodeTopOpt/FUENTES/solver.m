function [d_u] = solver(d_u,free_df,StifMat,fextlod,fix_df,fixnodes)

%condest(StifMat(free_df,free_df))
%d_u
if ~isempty(fix_df)
    d_u(fix_df)=fixnodes(:,3);
    d_u(free_df,1) = StifMat(free_df,free_df)\(fextlod(free_df) - StifMat(free_df,fix_df)*d_u(fix_df));
else
    d_u(free_df,1) = StifMat(free_df,free_df)\fextlod(free_df);
end




%fprintf(1,'COND NUMBER %d \n',condest(StifMat(free_df,free_df)))
end

