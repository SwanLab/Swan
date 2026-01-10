function [Fenv radii evecs chi2] = PlotEnvelopeFailure(x,y,z,LEGENDS,nfig)
if nargin == 0
    load('tmp.mat')
      
end
METHOD_PLOT =1 ;

   

[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( [x,y,z],'') ;

Fenv = [v(1) v(4) v(5)
    v(4) v(2) v(6)
    v(5) v(6) v(3)];
figure(nfig)
hold on
xlabel(LEGENDS{1})
ylabel(LEGENDS{2})
zlabel(LEGENDS{3})
%plot3( x, y, z, '.r' );
grid on
if METHOD_PLOT == 1
    nn = 50 ;
    [xx yy zz] = ellipsoid(center(1),center(2),center(3),radii(1),radii(2),radii(3),nn) ;
    
    % ROTATION
    X = xx(:);
    Y = yy(:);
    Z = zz(:) ;
    XMAT = evecs*[X';Y';Z'];
    X = XMAT(1,:) ;
    Y = XMAT(2,:) ;
    Z = XMAT(3,:) ;
    xx = reshape(X,nn+1,[]  );
    yy = reshape(Y,nn+1,[]  );
    zz = reshape(Z,nn+1,[]  );
    %
    
    mesh(xx,yy,zz)
  %  shading interp
   % colormap(gray);
    
    hidden off
else
    
    % draw data
    
    
    
    
    grid on
    
    
    hold on;
    
    %draw fit
    mind = min( [ x y z ] );
    maxd = max( [ x y z ] );
    nsteps = 50;
    step = ( maxd - mind ) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ),...
        linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
    
    Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
        2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
        2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
    p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
    set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
    
     
end
hold off;
view( -70, 40 );
axis vis3d  ;
camlight;
lighting phong;

% colormap hsv
% %alpha(.4)
