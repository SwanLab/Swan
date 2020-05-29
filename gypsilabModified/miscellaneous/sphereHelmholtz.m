function u = sphereHelmholtz(name,bndCond,rho,k,X) 
% Copyright (c) 20015-2017, Matthieu Aussal, Ecole Polytechnique       
% GNU General Public License v3.0. 
% Computation of analytic field for helmoltz spherical scattering

% Sign convention
if k <= 0
    sgnk = @(v) v;
else
    sgnk = @(v) conj(v);
end 
k = abs(k);

% Spherical data
[~,theta,r] = cart2sph(X(:,1),X(:,2),X(:,3));
theta = pi/2 - theta;  r = r'; 
u = zeros(1,length(theta));

% Boundary condition
if strcmp(bndCond,'dir')
    jsph = @(n,z) sqrt(pi./(2.*z)) .* besselj(n+0.5,z);
    hsph = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,2,z);
    alpha = @(n,z) - jsph(n,z)./hsph(n,z);    
elseif strcmp(bndCond,'neu')
    jsph = @(n,z) sqrt(pi./(2.*z)) .* besselj(n+0.5,z);
    hsph = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,2,z);
    djsph = @(n,z) (n*jsph(n-1,z) - (n+1)*jsph(n+1,z))/(2*n+1) ;    
    dhsph = @(n,z) -(1i./z.^2 - djsph(n,z)*hsph(n,z))/jsph(n,z);
    alpha = @(n,z) - djsph(n,z)./dhsph(n,z);
end

% Infinite spherical radiation
if strcmp(name,'inf')   
    n = 0;
    Pn = legendre(n,cos(theta));
    add = alpha(n,k*rho) .* Pn(1,:);
    while norm(add,'inf') > 1e-12
        u = u + add;
        n = n + 1;
        Pn = legendre(n,cos(theta));
        add = (-1)^n * (2*n+1) * alpha(n,k*rho) .* Pn(1,:);
    end
    u = (1i/k) .* u.';
    
% Finite spherical radiation
elseif strcmp(name,'bnd') || strcmp(name,'dom')
    n = 0;
    while abs(alpha(n,k*max(r))) > 1e-12
        Pn = legendre(n,cos(theta));
        u = u + (1i)^n * (2*n+1) .* ...
            (alpha(n,k*rho).*hsph(n,k*r)) .* Pn(1,:); %
        n = n+1;
    end
    u = u.';
    if strcmp(name,'dom')
        u(r<rho) = 0;
    end
    
else
    error('sphereHelmholtz.m : unavailable case')
end

% output
u = sgnk(u);
end