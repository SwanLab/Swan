function [zs] = besselJRobinzeros(c,k,N,freqCenter)

if ~exist('freqCenter','var')
    freqCenter = 0;
end

xmin = max(freqCenter - pi*N,0);
xmax = freqCenter + 2*pi*N; % Using zn ~ pi*n
if c == Inf
    % Returns the N first zeros of the function
    % Jk(x)
    zs = AllZeros(@(x)(besselj(k,x)),xmin,xmax,5*(N+freqCenter));
    
else
    % Returns the N first zeros of the function
    % c Jk(x) + xJk'(x)
    
    
    switch N
        case 0
            derBess = @(x)(-besselj(1,x));
        otherwise
            derBess = @(x)(0.5*(besselj(k-1,x)-besselj(k+1,x)));
    end
    
    zs = AllZeros(@(x)(c*besselj(k,x)+x*derBess(x)),xmin,xmax,5*(N+freqCenter));
    
    
    
end

[~,I] = sort(abs(zs-freqCenter));
zs = zs(I);

if zs(1)==0
    zs = zs(2:(N+1));
else
    zs = zs(1:N);
end

zs = zs(:);

end

function z=AllZeros(f,xmin,xmax,N)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros
if (nargin<4)
    N=100;
end
options=optimset('Display','off');
dx=(xmax-xmin)/N;
x2=xmin;
y2=f(x2);
z=[];
for i=1:N
    x1=x2;
    y1=y2;
    x2=xmin+i*dx;
    y2=f(x2);
    if (y1*y2<=0)                              % Rolle's theorem : one zeros (or more) present
        z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2),options)]; % Linear approximation to guess the initial value in the [x1,x2] range.
    end
end

end