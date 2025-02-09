function u = diskHelmholtz(name,bndCond,rho,k,X) 
% Copyright (c) 2018, M. Aussal &  M. Averseng, Ecole Polytechnique       
% GNU General Public License v3.0. 
% Computation of analytic field for helmoltz disk scattering

% Spherical data
[theta,r] = cart2pol(X(:,1),X(:,2));
theta = pi/2 - theta;  
u = zeros(length(theta),1);

% Boundary condition
if strcmp(bndCond,'dir')
    jcyl = @(n,z) besselj(n,z);
    hcyl = @(n,z) besselh(n,1,z);
    alpha = @(n,z) - jcyl(n,z)./hcyl(n,z);    
elseif strcmp(bndCond,'neu')
    jcyl = @(n,z) besselj(n,z);
    hcyl = @(n,z) besselh(n,1,z);
    djsph = @(n,z) 0.5*(jcyl(n-1,z)-jcyl(n+1,z)); 
    dhsph = @(n,z) 0.5*(hcyl(n-1,z)-hcyl(n+1,z));
    alpha = @(n,z) - djsph(n,z)./dhsph(n,z);
end

% Infinite spherical radiation
if strcmp(name,'inf')   
    n   = 0;
    Pn  = cos(n*theta);
    add = alpha(n,k*rho) .* Pn;
    while (norm(add,'inf') > 1e-12)
        u   = u + add;
        n   = n + 1;
        en  = 1+(n>0);
        Pn  = cos(n*theta);
        add = (-1)^n * en * alpha(n,k*rho) .* Pn;
    end

% Finite spherical radiation
elseif strcmp(name,'bnd') || strcmp(name,'dom')
    n = 0;
    while (abs(alpha(n,k*max(r))) > 1e-12)
        Pn = cos(n*theta);
        en = 1+(n>0);
        u  = u + en .* (-1i)^n .* ...
            (alpha(n,k*rho).*hcyl(n,k*r)) .* Pn; %
        n = n+1;
    end
    if strcmp(name,'dom')
        u(r<rho) = 0;
    end
    
else
    error('diskHelmholtz.m : unavailable case')
end
end
