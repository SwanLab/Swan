% test ellipsoid fit

% create the test data:
% radii
a = 1;
b = 2;
c = 3;
[ s, t ] = meshgrid( 0 : 0.1 : pi, 0 : 0.3 : pi );
x = a * cos(s) .* cos( t );
y = b * cos(s) .* sin( t );
z = c * sin(s);
% rotation
ang = pi/6;
xt = x * cos( ang ) - y * sin( ang );
yt = x * sin( ang ) + y * cos( ang );
% translation
shiftx = 4;
shifty = 5;
shiftz = 6;
x = xt + shiftx;
y = yt + shifty;
z = z  + shiftz;

% add noise:
noiseIntensity = 0.1; %
dx = randn( size( s ) ) * noiseIntensity;
dy = randn( size( s ) ) * noiseIntensity;
dz = randn( size( s ) ) * noiseIntensity;
x = x + dx;
y = y + dy;
z = z + dz;
x = x(:);
y = y(:);
z = z(:);
fprintf( 'Simulated average data deviation: %.5f\n', sqrt( sum( dx(:).^2 + dy(:).^2 + dz(:).^2 ) / size( x, 1 ) ) );

% do the fitting
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( [ x y z ], '' );
fprintf( 'Ellipsoid center: %.5g %.5g %.5g\n', center );
fprintf( 'Ellipsoid radii: %.5g %.5g %.5g\n', radii );
fprintf( 'Ellipsoid evecs:\n' );
fprintf( '%.5g %.5g %.5g\n%.5g %.5g %.5g\n%.5g %.5g %.5g\n', ...
    evecs(1), evecs(2), evecs(3), evecs(4), evecs(5), evecs(6), evecs(7), evecs(8), evecs(9) );
fprintf( 'Algebraic form:\n' );
fprintf( '%.5g ', v );
fprintf( '\nAverage deviation of the fit: %.5f\n', sqrt( chi2 / size( x, 1 ) ) );
fprintf( '\n' );

% draw data
figure,
plot3( x, y, z, '.r' );
hold on;

%draw fit
mind = min( [ x y z ] );
maxd = max( [ x y z ] );
nsteps = 50;
step = ( maxd - mind ) / nsteps;
[ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
hold off;
set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
view( -70, 40 );
axis vis3d equal;
camlight;
lighting phong;