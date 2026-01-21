function mesh = meshCrossLattice(L,t_frame,h1,h2,n1,n2,nFrame)
% L = 1;
% t_frame = 0.02;
% h1 = 1.38;
% h2 = 1.38;
% 
% n1 = 7;
% n2 = 7;

P = createPoints(L,t_frame,h1,h2);


%% ---------------------------------------------------------
%  BUILD 24 BLOCKS
% ---------------------------------------------------------
blocks = {};

%% 1) Diamond (4 blocks)
blocks{end+1} = add_block(P(20,:), P(17,:), P(18,:),P(19,:),n1,n1, 'Diam_NE');
% blocks{end+1} = add_block(P(25,:), P(21,:), P(18,:),P(22,:),n1,n2, 'Diam_NW');
% blocks{end+1} = add_block(P(23,:), P(25,:), P(22,:),P(19,:),n1,n2, 'Diam_SW');
% blocks{end+1} = add_block(P(20,:), P(24,:), P(25,:),P(23,:),n1,n2, 'Diam_SE');

%% 2) Bar 1 arms (4 arms × 2 blocks = 8 blocks)

blocks{end+1} = add_block(P(9,:), P(10,:),P(17,:),P(20,:),n1,n2, 'B1_NE_inner');
blocks{end+1} = add_block(P(19,:), P(18,:), P(13,:),P(14,:),n1,n2, 'B1_NE_outer');

% blocks{end+1} = add_block(P(22,:), P(18,:), P(13,:),P(7,:),n1,n2, 'B1_NW_inner');
% blocks{end+1} = add_block(P(19,:), P(22,:), P(7,:),P(14,:),n1,n2, 'B1_NW_outer');

%% 3) Bar 2 arms (4 arms × 2 blocks = 8 blocks)

blocks{end+1} = add_block(P(17,:), P(11,:), P(12,:),P(18,:),n2,n1, 'B2_NE_inner');
blocks{end+1} = add_block(P(16,:), P(20,:), P(19,:),P(15,:),n2,n1, 'B2_NE_outer');

% blocks{end+1} = add_block(P(23,:),P(19,:), P(15,:),P(8,:),n1,n2, 'B2_NW_inner');
% blocks{end+1} = add_block(P(16,:), P(20,:),P(23,:),P(8,:),n1,n2, 'B2_NW_outer');

%% 4) Frame (4 quadrants × 3 blocks = 12 blocks)
% Q1
blocks{end+1} = add_block(P(5,:), P(1,:), P(26,:),P(9,:),nFrame,n1, 'Frame_Q1_top');
blocks{end+1} = add_block(P(26,:), P(27,:) ,P(10,:),P(9,:),n1,nFrame, 'Frame_Q1_right');
blocks{end+1} = add_block(P(27,:) , P(2,:),P(6,:),P(10,:),n1,nFrame, 'Frame_Q1_inner');
%
% % Q2
blocks{end+1} = add_block(P(6,:),P(2,:), P(28,:),P(11,:),nFrame,n1, 'Frame_Q2_left');
blocks{end+1} = add_block(P(11,:), P(28,:), P(29,:),P(12,:),nFrame,n1, 'Frame_Q2_inner');
%
% % Q3
blocks{end+1} = add_block(P(12,:),P(29,:), P(3,:),P(7,:),nFrame,n1, 'Frame_Q3_bot');
blocks{end+1} = add_block(P(13,:), P(7,:), P(3,:),P(30,:),n1,nFrame, 'Frame_Q3_left');
blocks{end+1} = add_block(P(14,:), P(13,:), P(30,:),P(31,:) ,n1,nFrame, 'Frame_Q3_inner');

blocks{end+1} = add_block(P(8,:),P(14,:), P(31,:),P(4,:),n1,nFrame, 'Frame_Q4_bot');
blocks{end+1} = add_block(P(32,:),P(15,:), P(8,:),P(4,:),nFrame,n1, 'Frame_Q4_right');
blocks{end+1} = add_block(P(33,:),P(16,:), P(15,:),P(32,:),nFrame,n1, 'Frame_Q4_inner');
blocks{end+1} = add_block(P(1,:), P(5,:), P(16,:),P(33,:),nFrame,n1, 'Frame_Q4_up');

mesh = createMesh(blocks,L);

end

function [n, e] = mesh_quad(p1, p2, p3, p4, n1, n2)
    s = linspace(0,1,n1+1);
    t = linspace(0,1,n2+1);
    [S, T] = meshgrid(s, t);

    X = (1-S).*(1-T)*p1(1) + S.*(1-T)*p2(1) + S.*T*p3(1) + (1-S).*T*p4(1);
    Y = (1-S).*(1-T)*p1(2) + S.*(1-T)*p2(2) + S.*T*p3(2) + (1-S).*T*p4(2);
    
    % 1. Transpose first so the columns represent the grid rows
    Xt = X'; 
    Yt = Y';

    % 2. Flatten into columns
    n = [Xt(:), Yt(:)];
    % ------------------------------------------------------------------

    e = zeros(n1*n2, 4); % Pre-allocate for speed
    k = 1;
    for j = 1:n2
        for i = 1:n1
            idx = (j-1)*(n1+1)+i;
            % This connectivity creates a CCW traversal in (S,T) space:
            % (i,j) -> (i+1,j) -> (i+1,j+1) -> (i,j+1)
            e(k,:) = [idx, idx+1, idx+n1+2, idx+n1+1];
            k = k + 1;
        end
    end
end

function mesh = createMesh(blocks, L)
    coord = [];
    connec = [];
    
    node_offset = 0;
    
    % 1. Stack all nodes from all blocks blindly
    for i = 1:numel(blocks)
        p = blocks{i}.p;
        [n, e] = mesh_quad(p(1,:), p(2,:), p(3,:), p(4,:), ...
            blocks{i}.n1, blocks{i}.n2);
        
        coord = [coord; n];
        connec = [connec; e + node_offset];
        node_offset = node_offset + size(n,1);
    end
    
    % 2. Merge duplicate nodes STABLY
    % 'uniquetol' sorts by coordinate value. We need to unsort it to match
    % the original creation order (topology first, geometry second).
    tol = 1e-10;
    
    % Get Unique sorted nodes (C), index of first appearance (ia), and map (ic)
    [C, ia, ic] = uniquetol(coord, tol, 'ByRows', true, 'DataScale', 1);
    
    % Sort 'ia' to find the order in which unique nodes actually appeared
    [~, perm] = sort(ia);
    
    % Reorder coordinates to match appearance order
    coord_clean = C(perm, :);
    
    % Create an inverse permutation map to fix connectivity
    % We map the Sorted Index -> Stable Index
    inv_perm = zeros(size(perm));
    inv_perm(perm) = 1:length(perm);
    
    % Remap connectivity: 
    % Original connectivity pointed to 'coord' -> 'ic' points to 'C' ->
    % 'inv_perm' points to 'coord_clean'
    connec_clean = inv_perm(ic(connec));
    
    % Ensure connectivity is shaped correctly
    connec_clean = reshape(connec_clean, size(connec));
    
    s.connec = connec_clean;
    s.coord = coord_clean;
    
    % 3. Apply Boundary Shifts
    % Fixed: Use tolerance check instead of '==' for floating point numbers
    delta = 1e-9;
    xmax = L; xmin = -L; ymax = L; ymin = -L;
    geo_tol = 1e-8; 
    
    % Top Right Corner
    mask = abs(s.coord(:,1) - xmax) < geo_tol & abs(s.coord(:,2) - ymax) < geo_tol;
    s.coord(mask, :) = s.coord(mask, :) + [-delta, -delta];
    
    % Bottom Right Corner
    mask = abs(s.coord(:,1) - xmax) < geo_tol & abs(s.coord(:,2) - ymin) < geo_tol;
    s.coord(mask, :) = s.coord(mask, :) + [-delta, +delta];
    
    % Top Left Corner
    mask = abs(s.coord(:,1) - xmin) < geo_tol & abs(s.coord(:,2) - ymax) < geo_tol;
    s.coord(mask, :) = s.coord(mask, :) + [+delta, -delta];
    
    % Bottom Left Corner
    mask = abs(s.coord(:,1) - xmin) < geo_tol & abs(s.coord(:,2) - ymin) < geo_tol;
    s.coord(mask, :) = s.coord(mask, :) + [+delta, +delta];
    
    mesh = Mesh.create(s);
end

function blk = add_block(p1,p2,p3,p4,n1,n2,name)
 P = [p1; p2; p3; p4];

    % signed area (positive = CCW)
    A = 0.5 * sum( ...
        P(:,1) .* P([2 3 4 1],2) - ...
        P([2 3 4 1],1) .* P(:,2) );

    if A < 0
        % flip orientation
        P = P([1 4 3 2],:);
    end
blk.name = name;
blk.p    = [p1; p2; p3; p4];
blk.n1   = n1;
blk.n2   = n2;
end

function P = createPoints(L,t_frame,h1,h2)
%% outer corners
P(1,:) = [ L,  L];
P(2,:) = [-L,  L];
P(3,:) = [-L, -L];
P(4,:) = [ L, -L];

%% inner frame corners
P(5,:) = P(1,:) + [-t_frame, -t_frame];
P(6,:) = P(2,:) + [ t_frame, -t_frame];
P(7,:) = P(3,:) + [ t_frame,  t_frame];
P(8,:) = P(4,:) + [-t_frame,  t_frame];

%% bar offsets
a = (h1/2)/sin(pi/4);   % bar 1 offset
b = (h2/2)/sin(pi/4);   % bar 2 offset
c = max(a,b);

%% bar–frame intersection points
P(9,:)  = [0,P(5,2)] + [h1/2,0];
P(10,:) = [0,P(5,2)] - [h1/2,0];

P(11,:) = [P(6,1),0] + [0, h2/2];
P(12,:) = [P(6,1),0] - [0, h2/2];

P(13,:) = [0,P(7,2)] - [h1/2,0];
P(14,:) = [0,P(7,2)] + [h1/2,0];

P(15,:) = [P(6,2),0] - [0, h2/2];
P(16,:) = [P(6,2),0] + [0, h2/2];

%% diamond (intersection of bar strips)
P(17,:) = [-h1/2,  h2/2];   % N
P(18,:) = [-h1/2, -h2/2];   % W
P(19,:) = [h1/2,  -h2/2];  % S
P(20,:) = [h1/2,   h2/2];   % E

%% midpoints on diamond edges
P(21,:) = (P(17,:) + P(18,:))/2;
P(22,:) = (P(18,:) + P(19,:))/2;
P(23,:) = (P(19,:) + P(20,:))/2;
P(24,:) = (P(20,:) + P(17,:))/2;

%% center
P(25,:) = [0,0];
P(26,:) = P(9,:) + [ 0,t_frame] ;
P(27,:) = P(10,:) + [  0, t_frame] ;
P(28,:) = P(11,:) + [ -t_frame, 0];
P(29,:) = P(12,:) + [ -t_frame, 0] ;
P(30,:) = P(13,:) + [ 0,- t_frame];
P(31,:) = P(14,:) + [ 0,-t_frame] ;
P(32,:) = P(15,:) + [ t_frame,0];
P(33,:) = P(16,:) + [ t_frame,0] ;
% delta = (a+b)/2;
% P(26,:) = P(1,:) + [ 0,-delta] ;
% P(27,:) = P(1,:) + [  -delta,0] ;
% P(28,:) = P(2,:) + [  delta,0];
% P(29,:) = P(2,:) + [  0, -delta] ;
% P(30,:) = P(3,:) + [ 0, delta];
% P(31,:) = P(3,:) + [ delta,0] ;
% P(32,:) = P(4,:) + [ -delta,0];
% P(33,:) = P(4,:) + [ 0,delta] ;
end

function plotPoints(P)

figure; hold on; axis equal; grid on;

% scatter points
scatter(P(:,1), P(:,2), 50, 'filled');

% label each node with its ID
for i = 1:size(P,1)
    text(P(i,1), P(i,2), sprintf(' %d', i), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom');
end

xlabel('X');
ylabel('Y');
title('Node coordinates with node IDs');

hold off;
end

