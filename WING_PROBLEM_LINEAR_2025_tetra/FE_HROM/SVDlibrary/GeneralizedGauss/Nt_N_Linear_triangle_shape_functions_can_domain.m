function   [A, dA,DATAOUT ]=    Nt_N_Linear_triangle_shape_functions_can_domain(xLIM,x,DATA)
% Construct N^T*N, where N = [1-xi-eta, xi, eta]
 % ------------------------------------------------------------------------
% JAHO, 16-Jan-2022
% --------------------------____---------------------------------------------
if nargin == 0 
    load('tmp.mat')
end

  
dA = [] ;  
DATAOUT  = [] ; 
nfun = 3; 
N = cell(1,nfun) ; 
N{1} = 1- x(:,1)-x(:,2) ; 
N{2} = x(:,1) ; 
N{3} = x(:,2) ; 

A = zeros(size(x,1),9) ; 
iacum = 1; 
for ifun = 1:nfun 
    for jfun = 1:nfun
        A(:,iacum) = N{ifun}.*N{jfun} ; 
        iacum = iacum + 1; 
    end
end 



%  
% PLOTloc =0;
% polORDER = DATA.PORDER ; 
% 
% ndim = size(xLIM,1) ;
% if ndim ~=2
%     error('Incorrect option')
% end
% Mpoints =   size(x,1) ;
% nmonomials = (1+polORDER)^ndim ;
% 
% DATA = DefaultField(DATA,'EVALUATE_GRADIENT',0) ;
% 
% [Ax,dAx,DATAOUT_x ]= LagrangePolynomial1D(polORDER,xLIM(1,:),x(:,1)) ;
% [Ay,dAy,DATAOUT_y ]= LagrangePolynomial1D(polORDER,xLIM(2,:),x(:,2)) ;
% 
% DATAOUT.xNODES = DATAOUT_x.xnodes ;
% DATAOUT.yNODES = DATAOUT_y.xnodes ;
% 
% A = zeros(nmonomials,size(Ax,2)) ;
% % ------------------------------------------------------------
% if DATA.EVALUATE_GRADIENT == 1
%     dA_x = A;
%     dA_y = A ;
%     H_xx = A;     H_xy = A;     H_yy = A;
% else
%     dA_x = [] ;     dA_y = [] ;     H_xx = [];
%     H_xy = [];     H_yy = [];
% end
% % -----------------------------------------------------------
% COORnodesELEMENTS = zeros(ndim,nmonomials) ;
% % -----------------------------------------------------------
% iacum = 1;
% for idimx = 1:size(Ax,1)
%     for idimy = 1:size(Ay,1)
%         A(iacum,:) =  Ax(idimx,:).*Ay(idimy,:) ;
%         COORnodesELEMENTS(1,iacum) = [DATAOUT.xNODES(idimx)] ;
%         COORnodesELEMENTS(2,iacum) = [DATAOUT.yNODES(idimy)] ;
%         if   DATA.EVALUATE_GRADIENT == 1
%             dA_x(iacum,:) =  dAx(idimx,:).*Ay(idimy,:) ;
%             dA_y(iacum,:) =  Ax(idimx,:).*dAy(idimy,:) ;
%             H_xx(iacum,:) =  DATAOUT_x.ddf(idimx,:).*Ay(idimy,:) ;
%             H_yy(iacum,:) =   Ax(idimx,:).*DATAOUT_y.ddf(idimy,:) ;
%             H_xy(iacum,:) =  dAx(idimx,:).*dAy(idimy,:) ;
%         end
%         iacum = iacum + 1;
%     end
% end
% DATAOUT.COORnodesELEMENTS = COORnodesELEMENTS;
% if PLOTloc == 1    
%     nnn = sqrt(size(x,1)) ;
%     MESHX = reshape(x(:,1),nnn,[])' ;
%     MESHY = reshape(x(:,2),nnn,[])' ;
%     figure(1678)
%     hold on
%     title(['Lagrange polynomial, order = ',num2str(polORDER)])
%     for ifun = 1:size(A,1)
%         funLOC = A(ifun,:) ;
%         funLOC = reshape(funLOC,nnn,[]) ;        
%                 surf(MESHX,MESHY,funLOC);        
%     end    
%     if DATA.EVALUATE_GRADIENT == 1   
%         figure(124)        
%         hold on
%         title('Derivative')
%         xlabel('x')
%         ylabel('y')
%         zlabel('dA_x')
%         for ifun = 1:size(A,1)
%             funLOC = dA_x(ifun,:) ;
%             funLOC = reshape(funLOC,nnn,[]) ;  
%             surf(MESHX,MESHY,funLOC)            
%         end       
%         figure(125)        
%         hold on
%         title('Second Derivative (H_{xx})')
%         xlabel('x')
%         ylabel('y')
%         zlabel('H_{xx}')
%         for ifun = 1:size(A,1)
%             funLOC = H_xx(ifun,:) ;
%             funLOC = reshape(funLOC,nnn,[]) ; 
%             surf(MESHX,MESHY,funLOC)            
%         end       
%         figure(126)        
%         hold on
%         title('Second Derivative (H_{yy})')
%         xlabel('x')
%         ylabel('y')
%         zlabel('H_{yy}')
%         for ifun = 1:size(A,1)
%             funLOC = H_yy(ifun,:) ;
%             funLOC = reshape(funLOC,nnn,[]) ; 
%             surf(MESHX,MESHY,funLOC)
%         end
%         figure(127)        
%         hold on
%         title('Second Derivative (H_{xy})')
%         xlabel('x')
%         ylabel('y')
%         zlabel('H_{xy}')
%         for ifun = 1:size(A,1)
%             funLOC = H_xy(ifun,:) ;
%             funLOC = reshape(funLOC,nnn,[]) ; 
%             surf(MESHX,MESHY,funLOC)
%         end
%     end
% end
% A = A' ;
% dA_x = dA_x' ;
% dA_y = dA_y' ;
% dA  = {dA_x,dA_y} ; 
% DATAOUT.H_xx = H_xx' ;
% DATAOUT.H_yy = H_yy' ;
% DATAOUT.H_xy = H_xy' ;
% 
% end
