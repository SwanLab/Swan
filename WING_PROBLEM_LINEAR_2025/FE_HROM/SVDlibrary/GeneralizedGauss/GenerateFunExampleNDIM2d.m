function [Xf,W,mu,xMAT,MPOINTS,XY] = GenerateFunExampleNDIM2d(M,P,DATA,DATAFUN)
% Snapshot matrix Xf and integration weights for the
% TEsting function: page 12 of "A ‘best points’ interpolation method for efficient approximation
% of parametrized functions", by N.C. Nguyen et al, 2007
%---------------------------------------------------
% % G(x;mu) = (1-x)*cos(3*pi*mu(x+1))*e^((-1+x)*mu)
% Domain
if nargin == 0
    load('tmp.mat')
end
xLIM = DATAFUN.xLIM;
Dmu = DATAFUN.muP ;
DATA = DefaultField(DATA,'PARTITIONED_Xf',0);
DATA = DefaultField(DATA,'PLOT_FUNCTION_SNAP',0);
DATAFUN = DefaultField(DATAFUN,'REPEATMATRIX',0);

% Number of integration points (spatial grid)
% Location of integration points (and associated weights)
DATA = DefaultField(DATA,'INITIAL_DISCRETIZATION_TENSORPRODUCT',0);
%dbstop('21')


% ----------------------------------------------------------
[MPOINTS, x, W] = InitalDiscretization(DATA,M,xLIM) ; 
% ndim = 2;
%dbstop('38')
XY = x; 
[xxGRID, yyGRID] = meshgrid(x{1},x{2}) ;
xx = xxGRID(:);
% Coordinate y
yy = yyGRID(:);
% Therefore
xMAT = [xx yy] ;
MESH{1} = [xxGRID] ; 
MESH{2} = [yyGRID] ; 
 
% -----------------------------------------------





% Snapshot matrix


% Grid for the parametric space
mu = linspace(Dmu(1),Dmu(2),P) ;

if DATAFUN.REPEATMATRIX>0
    f = DATAFUN.REPEATMATRIX ;
else
    f = 1 ;
end

disp('Generating snapshot matrix ...')
if DATAFUN.TYPE(1) ==10
    Pdim = (P+1)^2 ;
else
    Pdim = P ; 
end 
disp(['Size = ',num2str(f*M(1)*M(2)*Pdim*length(DATAFUN.TYPE)*8e-6),' Mbytes'])


M = length(W);

%if DATA.PARTITIONED_Xf == 0
Xf = [] ;
% XfDER_x = [] ;
% XfDER_y = [] ;
x = xMAT(:,1) ;
y = xMAT(:,2) ;

% 
% DATA.DATAREMOVEPOINTS.HOLE.TYPE = 'POLYGONAL'; 
% DATA.DATAREMOVEPOINTS.HOLE.POINTS = [-0.6 -0.6
%                                      0.6 -0.6
%                                      0.6 0.6
%                                      -0.6 0.6]; 
DATA.DATAREMOVEPOINTS = DefaultField(DATA.DATAREMOVEPOINTS,'HOLE',[] ) ;

DATA = DefaultField(DATA,'IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA',[])  ; % For irregular quadrilaterals 


for idim = 1:length(DATAFUN.TYPE)
    % dbstop('33')
    
    %     XflocDERx = zeros(M,P) ;
    %     XflocDERy = zeros(M,P) ;
    
    if  DATAFUN.TYPE(idim) ==1
        Xfloc = SinCosExpFun(mu,x,y) ;
    elseif DATAFUN.TYPE(idim) ==10
        Xfloc = Poly2Dfun(P,x,y) ;
    elseif DATAFUN.TYPE(idim) ==11
      if ~isempty(DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA)
          % See LagrangePolynomial2D_irregular_aux.mlx
          Xfloc =    LagrangePolynomial2D_irregular(xLIM,P,[x,y],DATA) ; 
      else
       Xfloc =    LagrangePolynomial2D(xLIM,P,[x,y]) ;
      end
      elseif DATAFUN.TYPE(idim) ==12 
        Xfloc =    LagrangePolynomial2D_B_B(xLIM,P,[x,y],DATA) ; 
        
         elseif DATAFUN.TYPE(idim) ==13 
             error('Option not verified properly')
        Xfloc =    LagrangePolynomial2D_N_N(xLIM,P,[x,y],DATA) ; 
      
    else
        error('Option not implemented')
    end
    Xf = [Xf Xfloc] ;
    %     XfDER_x =  [XfDER_x XflocDERx] ;
    %     XfDER_y =  [XfDER_y XflocDERy] ;
end

if ~isempty(DATA.DATAREMOVEPOINTS.HOLE)
    switch  DATA.DATAREMOVEPOINTS.HOLE.TYPE
        case 'POLYGONAL'
            POINTS_pol = DATA.DATAREMOVEPOINTS.HOLE.POINTS ;
            POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
            INPol = inpolygon(x,y,POINTS_pol(:,1),POINTS_pol(:,2)) ; 
            Xf(INPol,:) = 0 ; 
        otherwise 
            error('Option not implemented')
    end
end

% XfDER = {XfDER_x,XfDER_y} ;


% if DATAFUN.REPEATMATRIX > 0
%
%     Xf = repmat(Xf,1,DATAFUN.REPEATMATRIX ) ;
%
%
% end

%PARTITION =[] ;
% else
%     error('This part of the code needs revision')
%
%     switch DATA.PartitionMethod
%         case 'ROWWISE'
%             nrows =M ;
%             nrowsB = floor(nrows/ndom) ;
%             PARTITION = cell(ndom,1) ;
%             iacum = 1 ;
%             for i=1:ndom
%                 PARTITION{i} =[iacum:iacum+nrowsB-1];
%                 iacum = iacum+nrowsB ;
%             end
%             if  PARTITION{end}(end) ~=nrows
%                 PARTITION{end} = PARTITION{end}(1): nrows ;
%             end
%         case 'GIVEN'
%             error('implement this option !!!!!')
%         otherwise
%             error('option not implemented yet')
%     end
%
%     xALL = x;
%     Xf = cell(ndom,1) ;
%     for idom = 1:ndom
%         ROWS = PARTITION{idom} ;
%         x = xALL(ROWS);
%         Xf{idom } = zeros(length(ROWS),P);
%         for iMU = 1:P
%             expPART = exp((-1+x)*mu(iMU)) ;
%             Xf{idom}(:,iMU) = ((1-x).*cos(3*pi*mu(iMU)*(x+1)).*expPART)' ;
%
%         end
%     end
%     x = xALL ;
% end

