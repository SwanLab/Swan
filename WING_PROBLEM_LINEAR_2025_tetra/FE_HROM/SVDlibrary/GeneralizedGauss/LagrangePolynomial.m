function   [A,fDER] =    LagrangePolynomial(xLIM,polORDER,x)
% Construct Lagrange Polynomials  in cartesian domains
% INPUTS  xLIM --- Limits cartesian domain; size(xLIM) = [ndim,2]
% polORDER : Order of the polynomial  (p=0,1,2,....)
% x: Points at which the function is to be evaluated size(xGAUSS) =
% [Mpoints,ndim] 
% f = Matrix of polynomials.  size(f) = [Mpoints,nmonomials]
%  where nmonomials = (1+p)^ndim
% fDER = cell array with the derivatives of the functions 
% ------------------------------------------------------------------------
% JAHO, 12-Oct-2021
% --------------------------____---------------------------------------------
if nargin == 0
    xLIM = [-1,1]; 
    polORDER =10; 
    x = linspace(xLIM(1),xLIM(2),100)' ; 
end
ndim = size(xLIM,1) ; 
Mpoints =   size(x,1) ; 
nmonomials = (1+polORDER)^ndim ; 
fDER = [] ; 


if ndim == 1
    xLIMloc = xLIM ;
    A = LagrangePolynomial1D(polORDER,xLIMloc,x(:))' ;
    
    
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


