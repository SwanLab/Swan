function e = boundary_nodes(t)
%BOUNDARY_NODES Find boundary nodes E of triangulation T
    
% code from Per-Olof Persson <persson@berkeley.edu>

    edges = [t(:,[1,2]);
             t(:,[2,3]);
             t(:,[3,1])];
    edges = sort(edges, 2);
    [~, ix, jx] = unique(edges, 'rows');
    vec = histc(jx, 1:max(jx));
    qx = find(vec == 1);
    e = edges(ix(qx), :);
    e = unique(e(:));
end
