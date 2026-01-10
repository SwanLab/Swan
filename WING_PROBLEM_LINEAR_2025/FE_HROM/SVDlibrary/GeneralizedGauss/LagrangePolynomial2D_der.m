function   [Ax,Ay] =    LagrangePolynomial2D_der(xLIM,polORDER,x,DATALAGRANGE)
% Construct Derivative Lagrange Polynomials  in cartesian domains
% INPUTS  xLIM --- Limits cartesian domain; size(xLIM) = [ndim,2]
% polORDER : Order of the polynomial  (p=0,1,2,....)
% x: Points at which the function is to be evaluated size(xGAUSS) =
% [Mpoints,ndim]
% Ax = der(f)/dx , Ay = def(f)/dy
% ------------------------------------------------------------------------
% JAHO, 12-Oct-2021
% --------------------------____---------------------------------------------
if nargin == 0
    load('tmp1.mat')
    %     polORDER = 1 ;
    %     x = [0,0] ;
end
ndim = size(xLIM,1) ;
if ndim ~=2
    error('Incorrect option')
end
Mpoints =   size(x,1) ;
nmonomials = (1+polORDER)^ndim ;


Ax = LagrangePolynomial1D(polORDER,xLIM(1,:),x(:,1))' ;
Ay = LagrangePolynomial1D(polORDER,xLIM(2,:),x(:,2))' ;

PLOT1D = 0 ;
if PLOT1D == 1
    figure(456)
    hold on
    xlabel('x')
    ylabel('A_x')
    grid on
    for  ifun = 1:size(Ax,1)
        plot(x(:,1),Ax(ifun,:),'*')
        
        
    end
    
    
    figure(457)
    hold on
    grid on
    xlabel('y')
    ylabel('A_y')
    for  ifun = 1:size(Ay,1 )
        
        plot(x(:,2),Ay(ifun,:),'*')
        
        
    end
    
end

A = zeros(nmonomials,size(Ax,2)) ;
iacum = 1;
for idimx = 1:size(Ax,1)
    for idimy = 1:size(Ay,1)
        A(iacum,:) =  Ax(idimx,:).*Ay(idimy,:) ;
        iacum = iacum + 1;
    end
    
end

PLOTloc =0;
if PLOTloc == 1
    
    nnn = sqrt(size(x,1)) ;
    MESHX = reshape(x(:,1),nnn,[]) ;
    MESHY = reshape(x(:,2),nnn,[]) ;
    
    for ifun = 1:size(A,1)
        funLOC = A(ifun,:) ;
        funLOC = reshape(funLOC,nnn,[]) ;
        
        figure(1)
        hold on
        surf(MESHX,MESHY,funLOC)
        
    end
    
end


A = A' ; 


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


