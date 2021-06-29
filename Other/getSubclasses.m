function tb = getSubclasses(rootclass,rootpath)
% GETSUBCLASSES Display all subclasses
%
%   GETSUBCLASSES(ROOTCLASS, [ROOTPATH])
%       Lists all subclasses of ROOTCLASS and their node dependency.
%       ROOTCLASS can be a string with the name of the class, an
%       object or a meta.class().
%
%       It looks for subclasses in the ROOTPATH folder and in all
%       its subfolders (at any depth). If ROOTPATH is not specified,
%       it is set to the folder where ROOTCLASS is located.
%       ROOTCLASS can also be a negative integer indicating how many
%       folders above that of ROOTCLASS to begin the search.
%
%   TB = GETSUBCLASSES...
%       TB is a table with
%           .names  -  name of the subclass
%           .from   -  current node
%           .to     -  node of the direct subclass
%
%       If the function is called with no output, it will plot
%       a graph() with the dependencies.
%
%
% Example:
%
%   which sde
%   C:\Program Files\MATLAB\R2016a\toolbox\finance\finsupport\@sde\sde.m  % sde constructor
%
%   getSubclasses('sde','C:\Program Files\MATLAB\R2016a\toolbox\finance\finsupport\');
%
% See also: META.CLASS, GRAPH

rootclass = validateRootclass(rootclass);
if isempty(rootclass) || ~isa(rootclass, 'meta.class')
    error('getSubclasses:unrecognizedClass','Unrecognized class.')
end

if nargin < 2
    rootpath = 0;
end
rootpath = validateRootpath(rootpath, rootclass);

garray = addNewTree([],rootclass,{rootclass.Name});
garray = recurseFolder(rootpath, garray);
tb     = getFromToNodes(garray{1});

if nargout == 0
    G = graph(tb.from, tb.to, [], unique(tb.names,'stable'));
    G.plot()
end
end

% retrieve files and subfolders recursively
function garray = recurseFolder(folder, garray)

[files, subfld] = getEligible(folder);
garray          = checkFiles(files,garray);

for ii = 1:numel(subfld)
    garray = recurseFolder(subfld{ii}, garray);
    %     cell2table(garray{1}.values,'VariableNames',garray{1}.keys)
end
end

% List files and folders to search down for subclasses
function [files,subfolders] = getEligible(folder)
list = dir(folder);
idir = [list.isdir];

files = {list(~idir).name};
files = regexp(files,'.+(?=\.m|\.p)','match','once');

subfolders = {list(idir).name};

% Drop . and .. (setdiff is slow)
nfld = nnz(idir);
idir = false(nfld,1);
for ii = 1:nfld
    idir(ii) = subfolders{ii}(1) == '.';
end
subfolders(idir) = [];

subfolders = fullfile(folder, subfolders);
end

% check if files are classes and build up parent relationship
function garray = checkFiles(files,garray)
for ii = 1:numel(files)
    mcls = meta.class.fromName(files{ii});
    if isempty(mcls) || isempty(mcls.SuperclassList) || isInTreeArray(garray, mcls.Name)
        % Do something with classes with no parents?
        continue
    end
    garray = addNewTree(garray, mcls);
    [garray,trash,hasMerged] = getSuperclassTree(mcls, garray, numel(garray), 1);
    if hasMerged
        garray(end) = [];
    end
end
end

% recursively travel up through all parents assigning node values
function [garray, current_node, doMerge] = getSuperclassTree(mcls,garray,index, current_node)
disp(mcls.Name)
% Add direct parents
[garray{index}, parent_list, current_node] = addParents(mcls, garray{index}, current_node);

% Check if any parent is mergeable
doMerge  = false;
nparen   = numel(parent_list);
iexclude = false(nparen,1);
for ii = 1:nparen
    p               = parent_list(ii);
    [doMerge, into] = isInTreeArray(garray(1:end-1), p.Name);
    if doMerge
        garray{into} = mergeTrees(garray{into}, garray{index}, p.Name);
        iexclude(ii) = true;
    end
end

% Exclude merged parents from recursion
parent_list(iexclude) = [];

% Recurse each parent up
for p = parent_list(:)'
    [garray, current_node, doMerge] = getSuperclassTree(p, garray, index, current_node);
end
end

function [tree, list, nodenum] = addParents(mcls, tree, nodenum)
[list, parent_names] = getParentList(mcls);
for ii = 1:numel(list)
    nodenum                   = nodenum + 1;
    [trash,grandparent_names] = getParentList(list(ii));
    tree(parent_names{ii})    = makeVal(grandparent_names,nodenum);
end
end
% get meta classes and names of parents
function [list,names] = getParentList(mcls)
list = mcls.SuperclassList;
if isempty(list)
    list  = [];
    names = mcls.Name;
    return
end

if nargout == 2
    names = arrayfun(@(x) x.Name,list,'un',0);
end
end

% merge tree into parent containing the root class of interest
function main = mergeTrees(main, subtree, mergeat)

num_children = subtree.Count - 1;
shift        = getNode(main, mergeat);
for k = main.keys
    val = main(k{1});
    if val.Node > shift
        val.Node   = val.Node + num_children;
        main(k{1}) = val;
    end
end

nodenum = getNode(subtree, mergeat);
for k = subtree.keys
    val = subtree(k{1});
    if val.Node < nodenum
        val.Node   = val.Node + shift;
        main(k{1}) = val;
    end
end
end

function node = getNode(g, classname)
val  = g(classname);
node = val.Node;
end

% checks if class is already part of any tree
function [bool,into] = isInTreeArray(array,classname)
numtrees = numel(array);
bool     = false;
into     = 0;
for t = 1:numtrees
    bool = bool | array{t}.isKey(classname);
    if bool
        into = t;
        return
    end
end
end

% initialize data structure that builds parent relationships
function array = addNewTree(array, class, parents)
if nargin < 3
    [trash, parents] = getParentList(class);
    array            = addNewTree(array, class, parents);
else
    array = [array, {containers.Map(class.Name, makeVal(parents,1))}];
end
end

function val = makeVal(parents, nodenum)
if ischar(parents)
    parents = {parents};
end
val = struct('Parent',{parents},'Node',nodenum);
end

% create list with from-to nodes for graph() input
function tb = getFromToNodes(g)
keys      = g.keys;
n         = numel(keys);
c         = 0;
[from,to] = deal(zeros(n,1));
names     = cell(n,1);
for k = keys
    value = g(k{1});
    for p = value.Parent(:)'
        c = c+1;
        try
            to(c) = getNode(g,p{1});
        catch
            to(c) = value.Node;
        end
            from(c)  = value.Node;
            names{c} = k{1};
    end
end
tb = table(names, from, to);
tb = unique(tb,'rows');
tb = sortrows(tb,'from');

% remove duplicates
[un,trash,subs]             = unique(tb.names);
irepeated                   = accumarray(subs,1) > 1;
irepeated                   = ismember(tb.names,un(irepeated));
icircular                   = tb.from == tb.to;
tb(irepeated & icircular,:) = [];

% TODO: avoid re-mapping node
[un,trash,from_new] = unique(tb.from);
[trash,pos]         = ismember(tb.to, tb.from);
to_new              = from_new(pos);
tb.from             = from_new;
tb.to               = to_new;
end

function rclass = validateRootclass(rclass)
if ischar(rclass)
    rclass = meta.class.fromName(rclass);
elseif isobject(rclass) && ~isa(rclass, 'meta.class')
    rclass = meta.class.fromName(class(rclass));
end
end

function rpath = validateRootpath(rpath, rclass)
if isnumeric(rpath) && rpath <= 0 && mod(rpath,1) == 0
    s = what(rclass.Name);
    for ii = 1:abs(rpath)
        s.path = fileparts(s.path);
    end
    rpath = s.path;
elseif ~ischar(rpath)
    error('getSubclasses:unrecognizedPath','Unrecognized path.')
end
end