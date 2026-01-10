function [wGAUSS,xGAUSS,errorGAUSS] = StandardGaussRule(DATA,EXACT_INTEGRAL,ndim,DATAFITTING)

if nargin == 0
    load('tmp1.mat')
    DATA.NUMBER_GAUSS_EACH_DIRECTION =[60,60] ; 
end

DATA = DefaultField(DATA,'NUMBER_GAUSS_ONE_DIRECTION',[]);
DATA = DefaultField(DATA,'NUMBER_GAUSS_EACH_DIRECTION',[]);

DETERMINE_GAUSS = 1;
DATA = DefaultField(DATA,'TYPEFUN',[]);


if isempty(DATA.NUMBER_GAUSS_ONE_DIRECTION) && ~isempty(DATA.TYPEFUN)
    NAMEFUNTYPE =  DATA.TYPEFUN  ;
    switch NAMEFUNTYPE
        case {'LagrangePolynomial3Dgen','LagrangePolynomial2Dgen','LagrangePolynomial1Dgen'}
            P = DATA.PORDER ;
            m_1d = ceil((P+1)/2) ;
        case 'LagrangePolynomial_ALLD_B_B'
            P = DATA.PORDER ;
            m_1d = ceil((2*P+1)/2) ;
            
            %         case 'Fourier1Dgen'
            %             m_1d = DATA.PORDER + 1;
            
        otherwise
            m_1d = DATA.PORDER_GAUSS ;
            % error('Specify the number of Gauss points (DATA.PORDER_GAUSS)')
    end
    m_1d = m_1d*ones(1,ndim) ;
elseif ~isempty(DATA.NUMBER_GAUSS_EACH_DIRECTION)
    m_1d = DATA.NUMBER_GAUSS_EACH_DIRECTION ;
    
elseif ~isempty(DATA.NUMBER_GAUSS_ONE_DIRECTION)
    
    m_1d = DATA.NUMBER_GAUSS_ONE_DIRECTION ;
    m_1d = m_1d*ones(1,ndim) ;
else
    DETERMINE_GAUSS = 0 ;
end

if DETERMINE_GAUSS==0
    wGAUSS= [ ]; xGAUSS = [] ; errorGAUSS = [] ;
else
    if ndim == 3
        [MPOINTS, Mglo, x, wGAUSS]  =TensorProd3Ddiscr(m1d,DATA.xLIM) ;
        
        [xx, yy,zz] = meshgrid(x{1},x{2},x{3}) ;
        xx = xx(:);
        yy = yy(:);
        zz = zz(:) ;
        xGAUSS = [xx yy zz] ;
    elseif ndim == 2
        [MPOINTS, Mglo, x, wGAUSS]  =TensorProd2Ddiscr(m_1d,DATA.xLIM) ;
        
        [xx, yy] = meshgrid(x{1},x{2}) ;
        xx = xx(:);
        yy = yy(:);
        
        xGAUSS = [xx yy] ;
        
    elseif ndim ==1
        switch  DATA.TYPEFUN
            case 'Fourier1Dgen'
                dX = DATA.xLIM(2)-DATA.xLIM(1) ;
                N = DATA.PORDER;
                wMID = dX/N ;
                nPOINTS = N + 1;
                wGAUSS = wMID*ones(nPOINTS,1) ;
                wGAUSS(1) = wGAUSS(1)/2 ; wGAUSS(end) = wGAUSS(end)/2;
                xGAUSS = linspace(DATA.xLIM(1),DATA.xLIM(2),nPOINTS)' ;
                
            otherwise
                [xGAUSS, wGAUSS] = GaussQuad(m_1d, DATA.xLIM(1), DATA.xLIM(2)) ;
                
        end
    else
        error('Option not implemented ')
        
    end
    [g]= Evaluate_Basis_Grad_Analytic(xGAUSS,[],DATA,0,DATAFITTING);
    
    
    APPR_INTEGRAL  = g*wGAUSS ;
    
    errorGAUSS = norm(EXACT_INTEGRAL-APPR_INTEGRAL)./norm(EXACT_INTEGRAL)*100 ;
    disp('-------------------------------------------------------------------------------------------')
    disp(['Integration Error (original function) using  STANDARD GAUSSIAN RULE with  ',num2str(length(wGAUSS)),' points = ',num2str(errorGAUSS),' %']) ;
    disp('-------------------------------------------------------------------------------------------')
end