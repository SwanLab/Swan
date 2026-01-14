function DATA = PointsIncludedECMcartesianGEN(DATA,x)


DATA = DefaultField(DATA,'PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM',[]) ;
DATA.ECM_POINTS_INCLUDE =[] ;
if ~isempty(DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM)
    %ndim = size(x,2) ;
    %   nperc_direct = ((DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM)^(1/ndim));
    nperc_include = 100-DATA.PERCENTAGE_POINTS_EXCLUDED_BOUNDARIES_DECM ;
    xLIM = DATA.xLIM ;
    xCENT =sum(xLIM,2)/2 ;
    xSPAN =  (xLIM(:,2)-xLIM(:,1))/2;
    
    xSPANincludeHALF = xSPAN*nperc_include/100;
    
    xINCLUDE_min = xCENT -xSPANincludeHALF ;
    xINCLUDE_max =  xCENT +xSPANincludeHALF ;
    xINC = [xINCLUDE_min,xINCLUDE_max] ;
    
    ndim = size(xLIM,1);
    IINN = ones(size(x,1),1) ;
    for idim = 1:ndim
        IINN = IINN.*(x(:,idim) <= xINC(idim,2)).*(x(:,idim) >= xINC(idim,1)) ;
    end
    
    
    IndInclude = find(IINN~=0) ;
    DATA.ECM_POINTS_INCLUDE= IndInclude ;
    
    
end

DATA = DefaultField(DATA,'ListGaussToExclude',[]) ; 

if ~isempty(DATA.ListGaussToExclude)
    DATA.ECM_POINTS_INCLUDE = setdiff(DATA.ECM_POINTS_INCLUDE,DATA.ListGaussToExclude);
end

