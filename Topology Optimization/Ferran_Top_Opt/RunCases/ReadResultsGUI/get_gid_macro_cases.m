function [numsym,corners_str,symcase,view_file,case_name] = get_gid_macro_cases(case_name)

%% Special cases (legacy compatibility)
switch case_name
    case 'RVE04N3'
        case_name = 'RVE_Square';
    case 'RVE04N3Fine'
        case_name = 'RVE_Square_Fine';
end

%% Load mesh
pointload = [];
sideload = [];
nodesolid = [];
dirichlet_data = [];
corners = [];
Group = [];
Initial_holes = [];
External_border_nodes = [];

% FICHERO CON LAS COORDENADAS Y CONECTIVIDADES
imesh = 1;
fname = [case_name,'_MESH_',num2str(imesh)];
eval(fname); 
[npt,ncol] = size(coord);
coordinates = coord(1:npt,2:ncol-1);

%% Number of symmetries
if ~isempty(strfind(case_name,'CantiliberbeamSym'))
    numsym = 0;
    view_file = 'CantiliberbeamSymFineFine.vv';
elseif ~isempty(strfind(case_name,'BridgeSimetric'))
    numsym = 1;
    view_file = 'BridgeSimetric.vv';
elseif ~isempty(strfind(case_name,'Gripping'))
    numsym = 0;
    view_file = 'GrippingFineFine.vv';
elseif ~isempty(strfind(case_name,'Square')) % also for RVE_Square
    numsym = 2;
    view_file = 'Square.vv';
elseif ~isempty(strfind(case_name,'Hexagonal'))
    numsym = 3;
    view_file = 'RVE_Hexagonal.vv';
elseif ~isempty(strfind(case_name,'Bicycle'))
    numsym = 0;
    view_file = 'Bicycle.vv';
elseif ~isempty(strfind(case_name,'TrencalosSupport'))
    numsym = 0;
    view_file = 'TrencalosSupport.vv';
elseif ~isempty(strfind(case_name,'CantileverExampleProblem'))
    numsym = 0;
    view_file = 'CantileverExampleProblem.vv';
else
    error('View not included yet (%s), please add it.',case_name);
end

%% Get corners
if numsym < 3
    maxcoord = max(coordinates);
    mincoord = min(coordinates);
    [~,corners(1)] = ismember(mincoord,coordinates,'rows');
    [~,corners(2)] = ismember([mincoord(1),maxcoord(2)],coordinates,'rows');
    [~,corners(3)] = ismember(maxcoord,coordinates,'rows');
    [~,corners(4)] = ismember([maxcoord(1),mincoord(2)],coordinates,'rows');
else
    [corners,~] = compute_periodic_boundary_nodes_hexagonal(coordinates);
end

if ~isempty(strfind(case_name,'BridgeSimetric')) % for these cases the symmetry line is on the right side
    corners = fliplr(corners);
end

% Correction for macro simplicity
nmax = 6; % maximum number of corners (from case RVE_Hexagonal)
corners(end+1:nmax) = 0;

corners_str = '';
for i = 1:length(corners)
    corners_str = [corners_str,' ',num2str(corners(i))];
end
corners_str = strtrim(corners_str);


%% Symmetry case
symcase = '1';
if ~isempty(strfind(case_name,'RVE'))
    symcase = '0';
end


end