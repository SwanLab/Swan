function z = genop(op,x,y)
%GENOP Generalized array operations.
%   GENOP(OP, X, Y) applies the function OP to the arguments X and Y where
%   singleton dimensions of X and Y have been expanded so that X and Y are
%   the same size, but this is done without actually copying any data.
%
%   OP must be a function handle to a function that computes an
%       element-by-element function of its two arguments.
%
%   X and Y can be any numeric arrays where non-singleton dimensions in one
%       must correspond to the same or unity size in the other.  In other
%       words, singleton dimensions in one can be expanded to the size of
%       the other, otherwise the size of the dimensions must match.
%
%   For example, to subtract the mean from each column, you could use
%
%       X2 = X - repmat(mean(X),size(X,1),1);
%
%   or, using GENOP,
%
%       X2 = genop(@minus,X,mean(X));
%
%   where the single row of mean(x) has been logically expanded to match
%   the number of rows in X, but without actually copying any data.
%
%   GENOP(OP) returns a function handle that can be used like above:
%
%       f = genop(@minus);
%       X2 = f(X,mean(X));

% written by Douglas M. Schwarz
% email: dmschwarz (at) urgrad (dot) rochester (dot) edu
% 13 March 2006

% This function was inspired by an idea by Urs Schwarz (no relation) and
% the idea for returning a function handle was shamelessly stolen from
% Duane Hanselman.

% Check inputs.
if ~(nargin == 1 || nargin == 3)
	error('genop:zeroInputs','1 or 3 arguments required.')
end
if ~isa(op,'function_handle')
	error('genop:incorrectOperator','Operator must be a function handle.')
end
if nargin == 1
	z = @(x,y) genop(op,x,y);
	return
end

% Compute sizes of x and y, possibly extended with ones so they match
% in length.
nd = max(ndims(x),ndims(y));
sx = size(x);
sx(end+1:nd) = 1;
sy = size(y);
sy(end+1:nd) = 1;
dz = sx ~= sy;
dims = find(dz);
num_dims = length(dims);

% Eliminate some simple cases.
if num_dims == 0 || numel(x) == 1 || numel(y) == 1
	z = op(x,y);
	return
end

% Check for dimensional compatibility of inputs, compute size and class of
% output array and allocate it.
if ~(all(sx(dz) == 1 | sy(dz) == 1))
	error('genop:argSizeError','Argument dimensions are not compatible.')
end
sz = max([sx;sy]);
z1 = op(x(1),y(1));
if islogical(z1)
	z = repmat(logical(0),sz);
else
	z = zeros(sz,class(z1));
end

% The most efficient way to compute the result seems to require that we
% loop through the unmatching dimensions (those where dz = 1), performing
% the operation and assigning to the appropriately indexed output.  Since
% we don't know in advance which or how many dimensions don't match we have
% to create the code as a string and then eval it.  To see how this works,
% uncomment the disp statement below to display the code before it is
% evaluated.  This could all be done with fixed code using subsref and
% subsasgn, but that way seems to be much slower.

% Compute code strings representing the subscripts of x, y and z.
xsub = subgen(sy ~= sz);
ysub = subgen(sx ~= sz);
zsub = subgen(dz);

% Generate the code.
indent = 2; % spaces per indent level
code_cells = cell(1,2*num_dims + 1);
for i = 1:num_dims
	code_cells{i} = sprintf('%*sfor i%d = 1:sz(%d)\n',indent*(i-1),'',...
		dims([i i]));
	code_cells{end-i+1} = sprintf('%*send\n',indent*(i-1),'');
end
code_cells{num_dims+1} = sprintf('%*sz(%s) = op(x(%s),y(%s));\n',...
	indent*num_dims,'',zsub,xsub,ysub);
code = [code_cells{:}];

% Evaluate the code.
% disp(code)
eval(code)


function sub = subgen(select_flag)
elements = {':,','i%d,'};
selected_elements = elements(select_flag + 1);
format_str = [selected_elements{:}];
sub = sprintf(format_str(1:end-1),find(select_flag));
