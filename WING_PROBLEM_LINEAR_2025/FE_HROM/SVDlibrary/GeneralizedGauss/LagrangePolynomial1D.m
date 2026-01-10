function [f,df,DATAOUT] = LagrangePolynomial1D(polORDER,xLIMloc,x)
% See LagrangePolynomial1D_aux.mlx

if nargin == 0
    load('tmp2.mat')
end
PLOT_LOCAL = 0;

nnodes =polORDER + 1;
xnodes = linspace(xLIMloc(1),xLIMloc(2),nnodes) ;
Mpoints =   size(x,1) ;

DATAOUT.xnodes = xnodes ;

f = ones(Mpoints,nnodes) ;

for i = 1:(polORDER+1)
    for j =1:(polORDER+1)
        if i ~= j
            f(:,i) = f(:,i).*(x-xnodes(j))/(xnodes(i)-xnodes(j)) ;
        end
    end
end


if  PLOT_LOCAL == 1
    %  figure(45)
  
    subplot(3,1,1)
      xlabel('x')
    ylabel('f')
    hold on
    grid on
    for iplot = 1:size(f,2)
        plot(x,f(:,iplot)) ;
    end
end

% Derivatives, see LagrangePolynomial1D_aux.mlx
df = zeros(Mpoints,nnodes) ;


% "DIRECT" VERSION
% ---------------
% for i =1: (polORDER+1)
%     for k=1:(polORDER+1)
%
%         % -----------------------------
%         if k ~=i
%              g = 1 ;
%             for j=1:(polORDER+1)
%                 if i~=j && j~=k
%                     g = g.*(x-xnodes(j))/(xnodes(i)-xnodes(j)) ;
%                 end
%             end
%              % -----------------------------
%         df(:,i) = df(:,i) +  g/(xnodes(i)-xnodes(k)) ;
%         end
%
%     end
% end

% Indirect way
for i =1: (polORDER+1)
    df(:,i) = LagrangianDerLOC(Mpoints,polORDER,xnodes,x,i) ;
end

if  PLOT_LOCAL == 1
    %   figure(45)
   
    subplot(3,1,2)
     xlabel('x')
    ylabel('df/dx')
    hold on
    grid on
    for iplot = 1:size(df,2)
        plot(x,df(:,iplot)) ;
    end
end




% Second derivative , see LagrangePolynomial1D_aux.mlx
%-----------------------------------------------------------
ddf = zeros(Mpoints,nnodes) ;

for i =1: (polORDER+1)
   for k = 1:(polORDER+1)
       if i~=k 
           R =  LagrangianDERlocSUB(Mpoints,polORDER,xnodes,x,i,k); 
           ddf(:,i) = ddf(:,i) + R/(xnodes(i)-xnodes(k)) ; 
       end
   end
end

if  PLOT_LOCAL == 1
    %   figure(45)
   
    subplot(3,1,3)
     xlabel('x')
    ylabel('d2f/dx2')
    hold on
    grid on
    for iplot = 1:size(ddf,2)
        plot(x,ddf(:,iplot)) ;
    end
end



f= f' ;
df = df' ;
ddf = ddf' ;

DATAOUT.ddf = ddf ;

end

function dfLOC = LagrangianDerLOC(Mpoints,polORDER,xnodes,x,i)

dfLOC = zeros(Mpoints,1) ;

for k=1:(polORDER+1)
    % -----------------------------
    if k ~=i
        g = 1 ;
        for j=1:(polORDER+1)
            if i~=j && j~=k
                g = g.*(x-xnodes(j))/(xnodes(i)-xnodes(j)) ;
            end
        end
        % -----------------------------
        dfLOC = dfLOC +  g/(xnodes(i)-xnodes(k)) ;
    end
end


end



function dfLOC = LagrangianDERlocSUB(Mpoints,polORDER,xnodes,x,i,k)

dfLOC = zeros(Mpoints,1) ;

for h=1:(polORDER+1)
    % -----------------------------
    if h ~=i && h ~=k
        g = 1 ;
        for j=1:(polORDER+1)
            if i~=j && j~=h && j~=k
                g = g.*(x-xnodes(j))/(xnodes(i)-xnodes(j)) ;
            end
        end
        % -----------------------------
        dfLOC = dfLOC +  g/(xnodes(i)-xnodes(h)) ;
    end
end


end