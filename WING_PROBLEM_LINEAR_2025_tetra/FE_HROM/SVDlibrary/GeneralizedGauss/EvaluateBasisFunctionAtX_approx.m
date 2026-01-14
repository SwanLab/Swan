function  [PHIk_y, dPHIk_y, POLYINFO]=     EvaluateBasisFunctionAtX_approx(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)

if nargin == 0
    load('tmp.mat')
end


DATA = DefaultField(DATA,'METHOD','LOCAL_FITTING');

 


switch DATA.METHOD
    
    case 'SVD_based_FITTING'
        
        POLYINFO = []; 
        
        npoints = size(xNEW,1) ;
        nfun = size(DATA.PHI,2) ;
        PHIk_y = zeros(npoints,nfun) ;
        dPHIk_y = cell(1,2) ;
        dPHIk_y{1} = zeros(npoints,nfun) ;
        dPHIk_y{2} = zeros(npoints,nfun) ;
        INDfitGLO = cell(npoints,1) ;
        COOR = DATA.COOR ;
        DATAS = DATA.SVD_INTERPOLATION;
        DATAS.DATAfitSVD.Vxpoly;
        
        for  iFUN = 1:length(DATAS.DATAfitSVD.Vxpoly)
            Vxpoly = DATAS.DATAfitSVD.Vxpoly{iFUN} ;
            Uypoly =DATAS.DATAfitSVD.Uypoly{iFUN}  ;
            Vxpoly_der = DATAS.DATAfitSVD.Vxpoly_der{iFUN}  ;
            Uypoly_der = DATAS.DATAfitSVD.Uypoly_der{iFUN}  ;
            S = DATAS.DATAfitSVD.S{iFUN} ;
            
            [f,df] = ...
                SplineSVDinterpolation_eval(xNEW,Vxpoly,Uypoly,Vxpoly_der,Uypoly_der,S,DATAS) ;
            PHIk_y(:,iFUN) = f ;
            dPHIk_y{1}(:,iFUN) = df{1} ;
            dPHIk_y{2}(:,iFUN) = df{2} ;
            
        end
        
        
        
    case 'LOCAL_FITTING'
        
        POLYINFO = [] ; 
        npoints = size(xNEW,1) ;
        nfun = size(DATA.PHI,2) ;
        PHIk_y = zeros(npoints,nfun) ;
        dPHIk_y = cell(1,2) ;
        dPHIk_y{1} = zeros(npoints,nfun) ;
        dPHIk_y{2} = zeros(npoints,nfun) ;
        INDfitGLO = cell(npoints,1) ;
        COOR = DATA.COOR ;
        % Local fitting through polynomial interpolation (2D)
        for ipoint = 1:npoints
            for ifun = 1:nfun
                F = DATA.PHI(:,ifun) ;
                
                xNEWloc = xNEW(ipoint,:) ;
                
                if ifun == 1
                    DATA.INDfit = [] ;
                else
                    DATA.INDfit = INDfitGLO{ipoint} ; % Index points for fitting
                end
                
                [f_approx,d_approx,INDfit] = EstimateFunAndDerivative2D(COOR,F,xNEWloc,DATA) ;
                if ifun == 1
                    INDfitGLO{ipoint} = INDfit ;
                end
                PHIk_y(ipoint,ifun) = f_approx ;
                dPHIk_y{1}(ipoint,ifun) = d_approx(1) ;
                dPHIk_y{2}(ipoint,ifun) = d_approx(2) ;
            end
        end
        
    case  {'FE_INTERPOLATION','FE_ONLY'}
        [PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionAtX_FEinterp(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)  ;
    otherwise
        
        error('Option not implemented')
        
end


