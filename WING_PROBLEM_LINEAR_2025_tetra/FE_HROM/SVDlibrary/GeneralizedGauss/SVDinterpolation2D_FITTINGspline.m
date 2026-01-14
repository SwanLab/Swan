function [DATAOUT] = SVDinterpolation2D_FITTINGspline(DATALOC,PHI,xx,yy)
% POLYNOMIAL FITTING OF THE SVD DECOMPOSITION OF EACH COLUMN OF PHI
% sHAHO, 3-APRIL-2020

if nargin == 0
    load('tmp1.mat')
end


DATALOC=  DefaultField(DATALOC,'TOL_SVD',1e-3) ;  % Tolerance for the SVD 
DATALOC=  DefaultField(DATALOC,'PLOT_MODES_PHI',0) ;
 
% Loop over columns of PHI  --- Transponse of G 

NFUN = size(PHI,2) ;
DATAOUT.Uypoly = cell(NFUN,1)  ;
DATAOUT.Vxpoly =cell(NFUN,1)  ;
DATAOUT.S= cell(NFUN,1) ;


Sall_numb = zeros(size(PHI,2),1) ;
Errorx_all = Sall_numb ;
Errory_all = Sall_numb ;
disp('Computing SVD fitting parameters ...')
for iCOL = 1:size(PHI,2)
    disp(['PHI(:,',num2str(iCOL),')']) ;
    PHIloc = PHI(:,iCOL) ;
    lx = length(xx);
    ly = length(yy) ;
    Fmat = reshape(PHIloc,ly,[]);
    DATALOC.NFUN = iCOL;
    DATALOC.PLOT_FUNCTIONS = DATALOC.PLOT_MODES_PHI ;
    [Uypoly,Vxpoly,Uypoly_der,Vxpoly_der,S,Errorx,Errory] = ...
        SplineSVDinterpolation_poly(xx,yy,Fmat,DATALOC) ;
    
    Sall_numb(iCOL) = length(S) ;
    
    DATAOUT.Uypoly{iCOL} = Uypoly;
    DATAOUT.Vxpoly{iCOL} = Vxpoly;
    DATAOUT.Uypoly_der{iCOL} = Uypoly_der;
    DATAOUT.Vxpoly_der{iCOL} = Vxpoly_der ;
    DATAOUT.S{iCOL} = S ;
    
    
    if ~isempty(Errorx)
        Errorx_all(iCOL) = Errorx;
        Errory_all(iCOL) = Errory;
        disp(['Fitting ERROR x =',num2str(Errorx)])
        disp(['Fitting ERROR x =',num2str(Errory)])
    end
end


figure(34)
subplot(2,1,1)
hold on
xlabel('Number of PHI column')
ylabel('Number of Singular Values ')
title('POLYNOMIAL FITTING SVD')
bar(Sall_numb)

subplot(2,1,2)
hold on
xlabel('Number of PHI column')
ylabel('Fitting error')
plot(Errorx_all)
plot(Errory_all)



 