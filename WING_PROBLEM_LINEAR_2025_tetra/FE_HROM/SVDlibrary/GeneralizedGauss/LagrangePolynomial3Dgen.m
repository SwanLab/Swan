function   [A, dA,DATAOUT ]=    LagrangePolynomial3Dgen(xLIM,x,DATALAGRANGE)
% Construct Lagrange Polynomials  in cartesian domains
% INPUTS  xLIM --- Limits cartesian domain; size(xLIM) = [ndim,2]
% polORDER : Order of the polynomial  (p=0,1,2,....)
% x: Points at which the function is to be evaluated size(xGAUSS) =
% [Mpoints,ndim]
% f = Matrix of polynomials.  size(f) = [Mpoints,nmonomials]
%  where nmonomials = (1+p)^ndim
% ------------------------------------------------------------------------
% JAHO, 12-Oct-2021/21-Oct-2021
% --------------------------____---------------------------------------------
if nargin == 0
    load('tmp2.mat')
   
end
ndim = size(xLIM,1) ;
if ndim ~=3
    error('Incorrect option')
end
polORDER = DATALAGRANGE.PORDER ; 
Mpoints =   size(x,1) ;
nmonomials = (1+polORDER)^ndim ;

% DATALAGRANGE = DefaultField(DATALAGRANGE,'COMPUTE_DERIVATIVE',1) ;
%
% if DATALAGRANGE.COMPUTE_DERIVATIVE == 1
%
% end

[Ax,dAx,DATAOUT_x  ]= LagrangePolynomial1D(polORDER,xLIM(1,:),x(:,1)) ;
[Ay,dAy,DATAOUT_y ]= LagrangePolynomial1D(polORDER,xLIM(2,:),x(:,2)) ;
[Az,dAz,DATAOUT_z  ]= LagrangePolynomial1D(polORDER,xLIM(3,:),x(:,3)) ;

DATAOUT.xNODES = DATAOUT_x.xnodes ;
DATAOUT.yNODES = DATAOUT_y.xnodes ;
DATAOUT.zNODES = DATAOUT_z.xnodes ;



 COORnodesELEMENTS = zeros(ndim,nmonomials) ;


A = zeros(nmonomials,size(Ax,2)) ;
 DATALAGRANGE = DefaultField(DATALAGRANGE,'EVALUATE_GRADIENT',0) ;
 if DATALAGRANGE.EVALUATE_GRADIENT == 1
     
     dA_x = A;      dA_y = A ;      dA_z = A ;
     H_xx = A;      H_xy = A;       H_xz = A;
     H_yy = A;   H_yz = A; 
     H_zz = A ; 
     
 else
      dA_x = [];      dA_y = [] ;      dA_z = [] ;
     H_xx = [];      H_xy = [];       H_xz = [];
     H_yy = [];   H_yz = []; 
     H_zz = []; 
 end




 iacum = 1;
 for idimx = 1:size(Ax,1)
     for idimy = 1:size(Ay,1)
         for  idimz = 1:size(Az,1)
             A(iacum,:) =  (Ax(idimx,:).*Ay(idimy,:)).*(Az(idimz,:)) ;
             
                 COORnodesELEMENTS(1,iacum) = [DATAOUT.xNODES(idimx)] ;
        COORnodesELEMENTS(2,iacum) = [DATAOUT.yNODES(idimy)] ;
        COORnodesELEMENTS(3,iacum) = [DATAOUT.zNODES(idimz)] ;
             
             if DATALAGRANGE.EVALUATE_GRADIENT == 1
                 
                 dA_x(iacum,:) =  (dAx(idimx,:).*Ay(idimy,:)).*Az(idimz,:) ;
                 
                 dA_y(iacum,:) =  (Ax(idimx,:).*dAy(idimy,:)).*Az(idimz,:) ;
                 
                 dA_z(iacum,:) =  (Ax(idimx,:).*Ay(idimy,:)).*dAz(idimz,:) ;
                 
                 H_xx(iacum,:) =  DATAOUT_x.ddf(idimx,:).*Ay(idimy,:).*Az(idimz,:) ;
                 H_yy(iacum,:) =   Ax(idimx,:).*DATAOUT_y.ddf(idimy,:).*Az(idimz,:) ;
                 H_zz(iacum,:) =   Ax(idimx,:).*Ay(idimy,:).*DATAOUT_z.ddf(idimz,:);
                 
                 H_xy(iacum,:) =  dAx(idimx,:).*dAy(idimy,:).*Az(idimz,:) ;  
                 H_xz(iacum,:) =  dAx(idimx,:).*Ay(idimy,:).*dAz(idimz,:) ;
                 
                 H_yz(iacum,:) =  (Ax(idimx,:).*dAy(idimy,:)).*dAz(idimz,:) ;
             end
             
             iacum = iacum + 1;
         end
     end
     
 end

% PLOTloc =0;
% if PLOTloc == 1
%     
%     nnn = sqrt(size(x,1)) ;
%     MESHX = reshape(x(:,1),nnn,[]) ;
%     MESHY = reshape(x(:,2),nnn,[]) ;
%     figure(1)
%     hold on
%     for ifun = 1:size(A,1)
%         funLOC = A(ifun,:) ;
%         funLOC = reshape(funLOC,nnn,[]) ;
%         
%         
%         surf(MESHX,MESHY,funLOC)
%         
%     end
%     
%     
%     figure(124)
%     hold on
%     xlabel('x')
%     ylabel('y')
%     zlabel('dA_x')
%     for ifun = 1:size(A,1)
%         funLOC = dA_x(ifun,:) ;
%         funLOC = reshape(funLOC,nnn,[]) ;
%         
%         
%         surf(MESHX,MESHY,funLOC)
%         
%     end
%     
% end


A = A' ;
dA_x = dA_x' ;
dA_y = dA_y' ;
dA_z = dA_z' ;

dA = {dA_x,dA_y,dA_z} ; 

DATAOUT.H_xx = H_xx' ;
DATAOUT.H_xy = H_xy' ;
DATAOUT.H_xz = H_xz' ;

DATAOUT.H_yy = H_yy' ;
DATAOUT.H_yz = H_yz' ;

DATAOUT.H_zz = H_zz' ;
DATAOUT.COORnodesELEMENTS = COORnodesELEMENTS;


%     PLOT_Fun = 1;
%     if PLOT_Fun ==1
%
%         figure(324)
%         hold on
%         grid on
%         xlabel('x')
%         ylabel('Polynomial')
%         for iii=1:size(A,2)
%             plot(x,A(:,iii))
%         end
%
%     end



end

%
% for idim = 1:size(xLIM,1)
%     % Loop over number of dimensions
%     xLIMloc = xLIM(idim,:) ;
%     nnodes =polORDER + 1;
%     xnodes = linspace(xLIMloc(1),xLIMloc(2),nnodes) ;
%
%     for ifun = 1:(polORDER+1)
%     end
%
% end


