function [esTheta, esPhi] = sphereMaxwell(radius, frequency, theta, phi, nMax)
% Compute the complex-value scattered electric far field of a perfectly
% conducting sphere using the mie series. Follows the treatment in
% Chapter 3 of 
%
% Ruck, et. al. "Radar Cross Section Handbook", Plenum Press, 1970.
%  
% The incident electric field is in the -z direction (theta = 0) and is
% theta-polarized. The time-harmonic convention exp(jwt) is assumed, and
% the Green's function is of the form exp(-jkr)/r.
% 
% Inputs:
%   radius: Radius of the sphere (meters)
%   frequency: Operating frequency (Hz)
%   theta: Scattered field theta angle (radians)
%   phi: Scattered field phi angle (radians)
%   nMax: Maximum mode for computing Bessel functions
% Outputs:
%   esTheta: Theta-polarized electric field at the given scattering angles
%   esPhi: Phi-polarized electric field at the given scattering angles
%
%   Output electric field values are normalized such that the square of the
%   magnitude is the radar cross section (RCS) in square meters.
%
%   Author: Walton C. Gibson, email: kalla@tripoint.org

% speed of light
c = 299792458.0;

% radian frequency
w = 2.0*pi*frequency;

% wavenumber
k = w/c;

% conversion factor between cartesian and spherical Bessel/Hankel function
s = sqrt(0.5*pi/(k*radius));

% mode numbers
mode = 1:nMax; 

% compute spherical bessel, hankel functions
[J(mode)] = besselj(mode + 1/2, k*radius); J = J*s;
[H(mode)] = besselh(mode + 1/2, 2, k*radius); H = H*s;
[J2(mode)] = besselj(mode + 1/2 - 1, k*radius); J2 = J2*s;
[H2(mode)] = besselh(mode + 1/2 - 1, 2, k*radius); H2 = H2*s;
 
% derivatives of spherical bessel and hankel functions
% Recurrence relationship, Abramowitz and Stegun Page 361
kaJ1P(mode) = (k*radius*J2 - mode .* J );
kaH1P(mode) = (k*radius*H2 - mode .* H );


% Ruck, et. al. (3.2-1)
An = -((i).^mode) .* ( J ./ H ) .* (2*mode + 1) ./ (mode.*(mode + 1));

% Ruck, et. al. (3.2-2), using derivatives of bessel functions 
Bn = ((i).^(mode+1)) .* (kaJ1P ./ kaH1P) .* (2*mode + 1) ./ (mode.*(mode + 1));
 
[esTheta esPhi] = mieScatteredField(An, Bn, theta, phi, frequency);
     
end


function [esTheta, esPhi] = mieScatteredField(An, Bn, theta, phi, frequency)

% Compute the complex-value scattered electric far field of a sphere
% using pre-computed coefficients An and Bn. Based on the treatment in:
%
% Ruck, et. al. "Radar Cross Section Handbook", Plenum Press, 1970.
%  
% The incident electric field is in the -z direction (theta = 0) and is
% theta-polarized. The time-harmonic convention exp(jwt) is assumed, and
% the Green's function is of the form exp(-jkr)/r.
% 
% Inputs:
%   An: Array of Mie solution constants
%   Bn: Array of Mie solution constants
%       (An and Bn should have the same length)
%   theta: Scattered field theta angle (radians)
%   phi: Scattered field phi angle (radians)
% Outputs:
%   esTheta: Theta-polarized electric field at the given scattering angles
%   esPhi: Phi-polarized electric field at the given scattering angles
%
%   Output electric field values are normalized such that the square of the
%   magnitude is the radar cross section (RCS) in square meters.
%
%   Author: Walton C. Gibson, email: kalla@tripoint.org

% speed of light
c = 299792458.0;
 
% wavenumber
k = 2.0*pi*frequency/c;

sinTheta = abs(sin(theta)); % note: theta only defined from from 0 to pi 
cosTheta = cos(theta); % ok for theta > pi 
 
% first two values of the Associated Legendre Polynomial
plm(1) = -sinTheta;
plm(2) = -3.0*sinTheta*cosTheta;

S1 = 0.0;
S2 = 0.0;

p = plm(1);

% compute coefficients for scattered electric far field
for iMode = 1:length(An)

        % derivative of associated Legendre Polynomial
    if abs(cosTheta) < 0.999999
        if iMode == 1
            dp = cosTheta*plm(1)/sqrt(1.0 - cosTheta*cosTheta);
        else
            dp = (iMode*cosTheta*plm(iMode) - (iMode + 1)*plm(iMode - 1))/sqrt(1.0 - cosTheta*cosTheta);
        end
    end
     
    if abs(sinTheta) > 1.0e-6
        term1 = An(iMode)*p/sinTheta;
        term2 = Bn(iMode)*p/sinTheta;
    end

    if cosTheta > 0.999999
        % Ruck, et. al. (3.1-12)
        val = ((i)^(iMode-1))*(iMode*(iMode+1)/2)*(An(iMode) - i*Bn(iMode));
        S1 = S1 + val;
        S2 = S2 + val;
    elseif cosTheta < -0.999999
        % Ruck, et. al. (3.1-14)
        val = ((-i)^(iMode-1))*(iMode*(iMode+1)/2)*(An(iMode) + i*Bn(iMode));
        S1 = S1 + val;
        S2 = S2 - val;
    else
        % Ruck, et. al. (3.1-6)
        S1 = S1 + ((i)^(iMode+1))*(term1 - i*Bn(iMode)*dp);
        % Ruck, et. al. (3.1-7)
        S2 = S2 + ((i)^(iMode+1))*(An(iMode)*dp - i*term2);
    end
    
    % recurrence relationship for next Associated Legendre Polynomial
    if iMode > 1
        plm(iMode + 1) = (2.0*iMode + 1)*cosTheta*plm(iMode)/iMode - (iMode + 1)*plm(iMode - 1)/iMode;
    end
    p = plm(iMode + 1);
end

% complex-value scattered electric far field, Ruck, et. al. (3.1-5)
esTheta = S1*cos(phi);
esPhi = -S2*sin(phi);

% normalize electric field so square of magnitude is RCS in square meters
esTheta = esTheta*sqrt(4.0*pi)/k;
esPhi = esPhi*sqrt(4.0*pi)/k;
end





