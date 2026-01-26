function [U,S,V,eSVD] = SVD_selective(SNAP,DATA,COLUMNS_RVEloc)


ROWDIV = size(SNAP,1);  ;
        COLDIV = cellfun(@length,COLUMNS_RVEloc) ;
        SNAP = mat2cell(SNAP,ROWDIV,COLDIV) ;
        epsilon = (DATA.REACTIONmodes_equal_DISPmodes==0) ;
        DATASVD.RELATIVE_SVD = 1 ;
        DATASVD.COMPLETE_SVD = 1;
        DATASVD.EPSILON_GLO = DATA.TOL_LOC;
        [U,S,V,ETIME,RankR,eSVD] = RSVDqpGEN(SNAP,epsilon,DATASVD) ;