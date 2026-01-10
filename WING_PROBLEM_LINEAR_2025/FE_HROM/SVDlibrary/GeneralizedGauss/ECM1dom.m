function [z,w,errorGLO,DATAOUT] = ECM1dom(Xf,W,DATA)
%
% Empirical cubature method
% -----------------------------
% Given Xf and W, EmpiricalCubatureMethod1dom returns a set of indices z
% and associated positive weights w so that the integration error of each
% column of Xf (as well as the integral of the domain)  is minimized: 
%     min(Xf'*W - Xf(z,:)'w  + (sum(W) -sum(w)))
%
% INPUTS
% ------
% Xf: M x P n matrix containing, columnwise, the integrand defined at the M
% integration points of the underlying grid (n is the number of entries of the  integrand, and P the number of
% samples)
% W: M x 1 vector containing the integration weights associated with the
% underlying grid 
%  Optional inputs 
% -------------------
%DATA.TOLsvdXf =1; % Truncation tolerance determining dominant modes of Xf (in %)
% of Xf . Default value = 1e-12
% DATA.npoints   = Number of integrations points to be selected (default value empty, which means that m = rank(J)) 
% DATA.TOL   = Tolerance for determining number of integration points (default value 1e-12)
% DATA. PLOT_ERROR, % If  = 1, the error versus number of points graph is plotted;

% 
% Joaquín A. Hernández, January 8-th - 2016
% See paper: Dimensional hyperreduction of nonlinear parameterized finite element models using optimized cubature
%  hernandez2016dimensional.pdf
%---------------------------------------
%
if nargin == 0
 
    load('J.mat')
    Xf = Jun' ; 
    W  =rand(size(Xf,1),1) ; 
    DATA = []  ;
    DATA.ImposeVolumeConstraint = 0; 
    
end

% Default inputs
DATA = DefaultField(DATA,'npoints',[]);
DATA = DefaultField(DATA,'TOLsvdXf',0);
DATA = DefaultField(DATA,'TOL',1e-12);
DATA = DefaultField(DATA,'PLOT_ERROR',1);
DATA = DefaultField(DATA,'PLOTsingVAL_Xf',1);
DATA = DefaultField(DATA,'PLOT_ERROR_INTEGRATION',1);
DATA = DefaultField(DATA,'bMIN_withONE',1); % New option
DATA = DefaultField(DATA,'PLOTsingVAL_Xf_All',DATA.PLOTsingVAL_Xf); % New option
DATA = DefaultField(DATA,'criterTOLecmACTUALERROR',1); % New option
DATA = DefaultField(DATA,'PlotPathAnalogy',0); % New option
DATA = DefaultField(DATA,'ImposeVolumeConstraint',1); % New option
DATA = DefaultField(DATA,'NO_ORTHOGONALIZATION_SNAPSHOTMATRIX',0); % New option
DATA = DefaultField(DATA,'IncludeWeights_Jmin',1); % New option


 
%DATA = DefaultField(DATA,'includeSingularValuesinJ',0); % New option
DATA = DefaultField(DATA,'maximumSIZExf',200); %

%%%
% Determining the size of the partitions
%dbstop('53')
% DATA = PartitionDecide(Xf,DATA) ; 

% 
% ndom= ceil(sizesnap/DATA.maximumSIZExf); 
% %dbstop('52')
% DATA.ndomainsSVD =ndom ; % 
% disp(['Number of domains= ',num2str(ndom)])
%

TIMETOTAL = tic ;

% if DATA.SVD_BEFORE_zeroaverage ==1
%     if  DATA.ndomainsSVD == 1
%     [Lambda,S,~,ERRORsvd] = SVDtruncated(Xf,DATA.TOLsvdXf,0);
% else
%     % dbstop('23')
%     DATA.PLOTsingVAL_Xf_All = 0 ;
%     DATA.TOLdomains = [] ;
%     
%     [Lambda,S,~,ERRORsvd] = SVD_partROW(Xf,DATA.ndomainsSVD,DATA,[],size(Xf,1)) ;
%     
% end
% end

% ------------------------
% Computation of J and b (Box 5.1, steps 3 to 5)
%dbstop('85')
% if DATA.NO_ORTHOGONALIZATION_SNAPSHOTMATRIX==1
%     
%     if DATA.IncludeWeights_Jmin == 1
%          J = bsxfun(@times,Xf,sqrt(W))' ;
%          b = J*sqrt(W) ; 
%     else
%         J = Xf';
%         b = J*W ;
%     end
%     Jnorm = sqrt(sum(J.*J,1)) ;
%     
% elseif     DATA.NO_ORTHOGONALIZATION_SNAPSHOTMATRIX==2
    
%     if DATA.IncludeWeights_Jmin == 1
%          J = bsxfun(@times,Xf,sqrt(W))' ;   
%             Jnorm2 = sqrt(sum(J.*J,2)) ;
%              J = bsxfun(@times,J,1./Jnorm2) ;
%             
%          b = J*sqrt(W) ; 
%     else
%         J = Xf';
%         b = J*W ;
%     end
%     Jnorm = sqrt(sum(J.*J,1)) ;
%     
% else
    [J,b,Jnorm,DATAOUT] = Jmin_bmin(Xf,W,DATA);
%end






% Computation of the reduced set of points z and associated weights w (Box 4.1)
%[z,w,errorGLO,kiter]= EmpiricalCubatureMethod(J,b,W,Jnorm,DATA,Xf) ; 
SingVal_F = [] ;  % NEW  FUNCTION (SVDlibrary, MArch 2020). 
DATAecm.TOL = DATA.TOL ; 
[z,w,errorGLO,DATAOUTecm]= EmpiricalCubatureMethod(J',SingVal_F,W,DATAecm)


 
TOTALTIME = toc(TIMETOTAL) ;

disp(['TOTALTIME = ',num2str(TOTALTIME)])

DATAOUT.J = J ; 
DATAOUT.b = b ; 
DATAOUT.kiter = DATAOUTecm.kiteraciones-1; 

% if  DATA.PlotPathAnalogy == 1
%     %dbstop('185')
%     PlotPathAnalogy(Xf,z,w,J,b,W,DATA)
% end

 
if DATA.PLOT_ERROR_INTEGRATION == 1
    %  dbstop('128')
    figure(701)
    hold on
    xlabel('Number of Points')
    ylabel('log(INTEG. ERROR)')
    plot(log10(errorGLO/100))
    
    
    figure(702)
    hold on
    xlabel('Number of Points')
    ylabel('(INTEG. ERROR) %')
    plot((errorGLO))
end

% [z,w,errorGLO,DATAOUT] = SecondECM(Xf(z,:),z,DATA,w) ; 

