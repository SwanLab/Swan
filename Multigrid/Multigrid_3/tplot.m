function tplot(p, t, u)
% TPLOT Plot triangular mesh P,T and, optionally, solution U
% code from Per-Olof Persson <persson@berkeley.edu>
    
    clf
    if nargin < 3
        patch('vertices',p, 'faces',t, 'facecol',[.8, 1, .8], 'edgecol','k');
    else
        patch('vertices',p, 'faces',t, 'facevertexcdata',u, ...
              'facecol','interp', 'edgecol','none');
        colorbar
    end
    set(gcf, 'renderer','zbuffer');
    axis equal
    drawnow
end
