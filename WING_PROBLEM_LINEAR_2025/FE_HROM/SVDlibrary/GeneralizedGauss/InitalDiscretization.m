function [MPOINTS, x, w ] = InitalDiscretization(DATA,M,xLIM)

if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'INITIAL_DISCRETIZATION_TENSORPRODUCT',0) ; 

if  DATA.INITIAL_DISCRETIZATION_TENSORPRODUCT ==0
    ndim = length(M) ;
    x= {} ;
    Mglo = 1;
    Wsing = 1 ;
    MPOINTS = [];
    for idim = 1:ndim
        np(idim) = 2*M(idim)+1 ;
        x{idim} = linspace(xLIM(idim,1),xLIM(idim,2),np(idim));
        % Integration weight (the same for all points)
        Wsing = Wsing*(x{idim}(3)-x{idim}(1));
        x{idim} = x{idim}(2:2:end-1)';
        Mglo = Mglo*length(x{idim}) ;
        MPOINTS(idim) = length(x{idim}) ;
    end
    w = Wsing*ones(Mglo,1) ;
else
    
    if length(M) == 2
        [MPOINTS, Mglo, x w]  =TensorProd2Ddiscr(M,xLIM) ; 
    elseif  length(M) == 3
       % dbstop('26')
         [MPOINTS, Mglo, x w]  =TensorProd3Ddiscr(M,xLIM) ; 
    elseif length(M) ==1 
        [x, w] = GaussQuad(M, xLIM(1), xLIM(2)) ;
        MPOINTS = M ; 
        x = {x} ;
    else 
        
        
        error('Option not implemented')
        
    end

    
end