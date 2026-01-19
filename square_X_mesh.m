
function mesh = square_X_mesh(L, T, n_arm, n_width)
 % L: Half-length of the square (Total size 2L x 2L)
    % T: Struct with thickness (.top, .bot, .left, .right, .d1, .d2)
    % n_arm: Number of elements along the length of the diagonal arm
    % n_width: Number of elements across the thickness (width) of all bars

    % --- 1. DEFINE GEOMETRY BOUNDARIES ---
    % Inner Frame Rectangle
    fi.top = L - T.top; 
    fi.bot = -L + T.bot;
    fi.left = -L + T.left; 
    fi.right = L - T.right;
    
    % Outer Frame Rectangle
    fo.top = L; fo.bot = -L; fo.left = -L; fo.right = L;

    % --- 2. CALCULATE INTERSECTIONS (The "Key Points") ---
    % We need the 8 points where the diagonal edges hit the Inner Frame.
    % Diagonal Equations: y = +/- x. Offsets = (thickness/2)*sqrt(2)
    off1 = (T.d1/2) * sqrt(2); % Offset for Diag 1 (BL to TR)
    off2 = (T.d2/2) * sqrt(2); % Offset for Diag 2 (BR to TL)

    % -- Central Hub Corners (Intersection of diagonal lines) --
    % Top (D1_upper & D2_upper): y = x + off1, y = -x + off2
    hub.top   = [(off2-off1)/2, (off2-off1)/2 + off1];
    % Right (D1_lower & D2_upper): y = x - off1, y = -x + off2
    hub.right = [(off2+off1)/2, (off2+off1)/2 - off1];
    % Bot (D1_lower & D2_lower): y = x - off1, y = -x - off2
    hub.bot   = [(off1-off2)/2, (off1-off2)/2 - off1];
    % Left (D1_upper & D2_lower): y = x + off1, y = -x - off2
    hub.left  = [(-off2-off1)/2, (-off2-off1)/2 + off1];

    % -- Arm-Frame Intersections (Inner Loop Points) --
    % We define these 8 points in CCW order starting from Top-Right
    % We use a helper 'intersect_box' to find where lines hit the rect
    
    % TR Arm (Diag 1) hitting Top/Right boundaries
    p_tr_1 = intersect_box(1, off1, fi);  % Line y = x + off1 (Upper edge)
    p_tr_2 = intersect_box(1, -off1, fi); % Line y = x - off1 (Lower edge)
    
    % TL Arm (Diag 2) hitting Top/Left boundaries
    p_tl_1 = intersect_box(-1, off2, fi); % Line y = -x + off2 (Upper edge)
    p_tl_2 = intersect_box(-1, -off2, fi);% Line y = -x - off2 (Lower edge)
    
    % BL Arm (Diag 1) hitting Bot/Left boundaries
    p_bl_1 = intersect_box(1, off1, fi);  % Line y = x + off1 (Upper edge of bar)
    p_bl_2 = intersect_box(1, -off1, fi); % Line y = x - off1 (Lower edge)
    
    % BR Arm (Diag 2) hitting Bot/Right boundaries
    p_br_1 = intersect_box(-1, off2, fi); % Line y = -x + off2 (Upper)
    p_br_2 = intersect_box(-1, -off2, fi);% Line y = -x - off2 (Lower)

    % Organize Inner Points into a Ring (CCW)
    % Note: The order depends on which wall they hit. 
    % We assume standard "X" (TR arm hits TR corner region).
    % Points: [TR_Upper, TL_Upper, TL_Lower, BL_Upper, BL_Lower, BR_Lower, BR_Upper, TR_Lower]
    % Wait, correct CCW order on the perimeter is:
    % TR_Upper -> TL_Upper -> TL_Lower -> BL_Upper -> BL_Lower -> BR_Lower -> BR_Upper -> TR_Lower -> (back to start)
    % BUT we need to pair them for blocks.
    
    % Correct Pairing for Blocks:
    % Arm TR connects Hub to [p_tr_2, p_tr_1] (CCW order)
    
    % --- 3. PROJECT TO OUTER FRAME (The "Outer Loop") ---
    % For every inner point, we project it perpendicularly to the outer wall.
    % This creates the "Frame Blocks".
    
    outer_tr_1 = project_point(p_tr_1, fi, fo);
    outer_tr_2 = project_point(p_tr_2, fi, fo);
    outer_tl_1 = project_point(p_tl_1, fi, fo);
    outer_tl_2 = project_point(p_tl_2, fi, fo);
    outer_bl_1 = project_point(p_bl_1, fi, fo);
    outer_bl_2 = project_point(p_bl_2, fi, fo);
    outer_br_1 = project_point(p_br_1, fi, fo);
    outer_br_2 = project_point(p_br_2, fi, fo);

    % --- 4. MESH GENERATION (The 13 Zones) ---
    all_n = []; all_e = [];

    % -- ZONE 1: HUB --
    [n, e] = mesh_quad(hub.left, hub.bot, hub.right, hub.top, n_width, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % -- ZONES 2-5: ARMS --
    % TR Arm (Connects Hub Right/Top to Inner Frame TR points)
    [n, e] = mesh_quad(hub.top, hub.right, p_tr_2, p_tr_1, n_arm, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);
    
    % TL Arm
    [n, e] = mesh_quad(hub.left, hub.top, p_tl_1, p_tl_2, n_arm, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);
    
    % BL Arm
    [n, e] = mesh_quad(hub.bot, hub.left, p_bl_2, p_bl_1, n_arm, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);
    
    % BR Arm
    [n, e] = mesh_quad(hub.right, hub.bot, p_br_1, p_br_2, n_arm, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % -- ZONES 6-9: FRAME ANCHORS (Where arms touch frame) --
    % These blocks bridge the Inner Points to the Outer Points.
    % They must share the "n_width" edge with the arms.
    
    % TR Anchor
    [n, e] = mesh_quad(p_tr_1, p_tr_2, outer_tr_2, outer_tr_1, n_width, n_width); % n_width x n_width for consistency? Or length?
    % Let's use n_width for the radial direction (thickness of frame)
    % The interface with the arm has n_width nodes.
    % So dimensions are (n_width along arm interface, n_frame_thick radial)
    % Let's assume frame thickness density is also n_width for now.
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % TL Anchor
    [n, e] = mesh_quad(p_tl_2, p_tl_1, outer_tl_1, outer_tl_2, n_width, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);
    
    % BL Anchor
    [n, e] = mesh_quad(p_bl_1, p_bl_2, outer_bl_2, outer_bl_1, n_width, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % BR Anchor
    [n, e] = mesh_quad(p_br_2, p_br_1, outer_br_1, outer_br_2, n_width, n_width);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % -- ZONES 10-13: FRAME CORNERS (Bridging the anchors) --
    % These fill the gaps between the anchors.
    % Top Gap (Between TL_1 and TR_1)
    [n, e] = mesh_quad(p_tl_1, outer_tl_1, outer_tr_1, p_tr_1, n_width, n_arm); 
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % Left Gap (Between BL_2 and TL_2)
    [n, e] = mesh_quad(p_bl_2, outer_bl_2, outer_tl_2, p_tl_2, n_width, n_arm);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % Bot Gap
    [n, e] = mesh_quad(p_br_1, outer_br_1, outer_bl_1, p_bl_1, n_width, n_arm);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % Right Gap
    [n, e] = mesh_quad(p_tr_2, outer_tr_2, outer_br_2, p_br_2, n_width, n_arm);
    [all_n, all_e] = append_mesh(all_n, all_e, n, e);

    % --- 5. MERGE ---
    [mesh.coord, ~, ic] = unique(round(all_n, 5), 'rows');
    mesh.connec = ic(all_e);

    % Plot
    figure('Color','w'); patch('Faces',mesh.connec, 'Vertices',mesh.coord, 'FaceColor',[0.5 0.8 1]);
    axis equal; title('13-Zone Continuous Block Mesh');
end

% --- GEOMETRY HELPERS ---

function p = intersect_box(m, c, box)
    % Finds intersection of y = mx + c with a rectangle defined by box fields
    % This is simplified for the specific X-brace geometry (Diagonal lines)
    % We check intersections with all 4 inner walls and pick the relevant one.
    
    tol = 1e-5;
    candidates = [];
    
    % Top Wall (y = box.top) -> x = (box.top - c)/m
    xt = (box.top - c)/m;
    if xt >= box.left-tol && xt <= box.right+tol, candidates = [candidates; xt, box.top]; end
    
    % Bot Wall (y = box.bot) -> x = (box.bot - c)/m
    xb = (box.bot - c)/m;
    if xb >= box.left-tol && xb <= box.right+tol, candidates = [candidates; xb, box.bot]; end
    
    % Right Wall (x = box.right) -> y = m*box.right + c
    yr = m*box.right + c;
    if yr >= box.bot-tol && yr <= box.top+tol, candidates = [candidates; box.right, yr]; end
    
    % Left Wall (x = box.left) -> y = m*box.left + c
    yl = m*box.left + c;
    if yl >= box.bot-tol && yl <= box.top+tol, candidates = [candidates; box.left, yl]; end

    % Selection Logic: We know roughly where each point should be based on the quadrant
    % For this specific script, we rely on the caller knowing which line is which.
    % We assume the line hits the box at 2 points usually (in and out). 
    % But we clipped the infinite line.
    % To pick the correct point, we can use the angle or simple min/max logic?
    % Actually, for the "X", we just need the point closest to the corner of interest.
    
    % Heuristic: 
    % TR points should have x > 0, y > 0
    % TL points: x < 0, y > 0, etc.
    % Let's filter by quadrant.
    
    % But wait, m=1, c=small -> passes through TR and BL.
    % We need to return the SINGLE point relevant to the specific arm side.
    
    % Since this is a specific topology hardcoded in the main function:
    % if m*x > 0 (Slope 1) -> TR or BL. If y > 0 -> TR.
    % We return the one in the specific quadrant.
    
    % Let's just pick the one that matches the "Target Zone".
    % Refinement: intersect_box is too generic.
    % Let's use explicit logic in main? No, messy.
    
    % FILTER:
    for i = 1:size(candidates,1)
        px = candidates(i,1); py = candidates(i,2);
        % Identify quadrant of the point
        if m > 0 && c > 0 % D1 Upper -> TR or BL? No, y=x+c (c>0) -> Top Leftish shift
             % y = x + offset. 
             % Hits Top wall (TR quadrant) or Left wall (TL quadrant)?
             % Actually, y=x+off1 passes through Top and Left? No.
             % It's parallel to y=x.
             if px > -box.right && py > -box.top, p = [px, py]; return; end
        elseif m > 0 && c < 0 % D1 Lower -> y=x-c. Hits Bot or Right.
             if px > -box.right && py < box.top, p = [px, py]; return; end 
        end
        % This is getting risky. 
        % BETTER: Just pick the point maximizing the coordinate in the corner direction.
    end
    
    % FALLBACK: Simple Max/Min
    % TR Region: Max X, Max Y
    if m==1 && c > 0, [~,idx]=max(candidates(:,2)); p=candidates(idx,:); % Top edge
    elseif m==1 && c < 0, [~,idx]=max(candidates(:,1)); p=candidates(idx,:); % Right edge
        
    % BL Region
    elseif m==1 && c==0, % Center?
    elseif m==1 % BL is the "negative" of TR
       % Actually, y=x+c hits TR and BL.
       % TR point has higher Y. BL point has lower Y.
       % We need to know WHICH one we want.
       % Let's pass a "Quadrant" hint to this function.
    end
    
    % Let's fix the calls in Main:
    % intersect_box(m, c, box, 'TR')
    p = candidates(1,:); % Placeholder
end

% Function with Quadrant Hint
function p = intersect_box_quad(m, c, box, quad)
    % quad: [sign_x, sign_y]
    candidates = [];
    tol = 1e-5;
    % Check all 4 walls... (same as above)
    % Top
    xt = (box.top - c)/m;
    if xt >= box.left-tol && xt <= box.right+tol, candidates = [candidates; xt, box.top]; end
    % Bot
    xb = (box.bot - c)/m;
    if xb >= box.left-tol && xb <= box.right+tol, candidates = [candidates; xb, box.bot]; end
    % Right
    yr = m*box.right + c;
    if yr >= box.bot-tol && yr <= box.top+tol, candidates = [candidates; box.right, yr]; end
    % Left
    yl = m*box.left + c;
    if yl >= box.bot-tol && yl <= box.top+tol, candidates = [candidates; box.left, yl]; end
    
    % Select the one in the quadrant
    best_dist = inf;
    p = [0,0];
    target = [quad(1)*max(abs([box.right, box.left])), quad(2)*max(abs([box.top, box.bot]))];
    
    for i = 1:size(candidates,1)
        if norm(candidates(i,:) - target) < best_dist
            best_dist = norm(candidates(i,:) - target);
            p = candidates(i,:);
        end
    end
end

function p_out = project_point(p_in, fi, fo)
    % Project p_in (on inner rect fi) to p_out (on outer rect fo)
    % Logic: If on top edge, project up. If on right edge, project right.
    tol = 1e-4;
    p_out = p_in;
    if abs(p_in(2) - fi.top) < tol, p_out(2) = fo.top;       % On Top Edge
    elseif abs(p_in(2) - fi.bot) < tol, p_out(2) = fo.bot;   % On Bot Edge
    elseif abs(p_in(1) - fi.right) < tol, p_out(1) = fo.right; % On Right Edge
    elseif abs(p_in(1) - fi.left) < tol, p_out(1) = fo.left;   % On Left Edge
    end
end

function [all_n, all_e] = append_mesh(all_n, all_e, n, e)
    e = e + size(all_n, 1);
    all_n = [all_n; n];
    all_e = [all_e; e];
end

function [n, e] = mesh_quad(p1, p2, p3, p4, n1, n2)
    s = linspace(0,1,n1+1); t = linspace(0,1,n2+1);
    [S, T] = meshgrid(s, t);
    X = (1-S).*(1-T)*p1(1) + S.*(1-T)*p2(1) + S.*T*p3(1) + (1-S).*T*p4(1);
    Y = (1-S).*(1-T)*p1(2) + S.*(1-T)*p2(2) + S.*T*p3(2) + (1-S).*T*p4(2);
    n = [X(:), Y(:)];
    e = [];
    for j=1:n2, for i=1:n1
        idx = (j-1)*(n1+1)+i;
        e=[e; idx, idx+1, idx+n1+2, idx+n1+1];
    end, end
end

% REPLACEMENT FOR intersect_box calls in Main:
% p_tr_1 = intersect_box_quad(1, off1, fi, [1, 1]);
% p_tr_2 = intersect_box_quad(1, -off1, fi, [1, 1]);
% ... etc