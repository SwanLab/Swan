function T = combinations(varargin)
%COMBINATIONS Generate all element combinations of arrays
%
%   T = COMBINATIONS(A1, A2, ..., An) returns all element combinations.
%   Each row in output table T is a combination with the first element
%   coming from A1, the second element coming from A2, and so on. The
%   number of rows in T is the product of the number of elements in each
%   input, and the number of variables in T is the number of inputs.
% 
%   The inputs can be vectors, matrices, or multidimensional arrays. If an
%   input Ak is not a vector, then COMBINATIONS treats it as a single
%   column vector, Ak(:).
%
%   Examples:
%      % All element combinations of two vectors
%      T = combinations([2 1], [5 4 6])
%
%      % All element combinations of inputs of mixed datatypes
%      A1 = [0 1; 2 3];
%      A2 = ["a" "b"];
%      A3 = categorical(["x" "y"]);
%      T = combinations(A1, A2, A3)
%
%   See also REPELEM, REPMAT, PERMS, NCHOOSEK, NDGRID, MESHGRID

%   Copyright 2022-2023 The MathWorks, Inc.

narginchk(1,Inf);

numVars = nargin;
varNumel = zeros(1,numVars);
input = varargin;
for k = 1:numVars
    varNumel(k) = numel(input{k});
end

m = prod(varNumel);

% Gracefully error like in perms()
[~,maxsize] = computer;
if m*numVars > maxsize
    error(message('MATLAB:pmaxsize'))
end

vars = cell(1,numVars);
varNames = cell(1,numVars);
noName = false(1,numVars);
for k = 1:numVars
    try
        var = input{k}(:);
    catch
        error(message('MATLAB:combinations:InvalidType',class(input{k})));
    end
    r = repelem(var,prod(varNumel(k+1:end)),1);
    vars{:,k} = repmat(r,prod(varNumel(1:k-1)),1);
    varNames{k} = inputname(k);
    noName(k) = isempty(varNames{k});
end

% Construct output table containing the combinations 
% Use default table variable names for inputs with no name
numNoName = sum(noName);
if numNoName > 0
    vdim = matlab.internal.tabular.private.varNamesDim;
    varNames(noName) = vdim.dfltLabels(find(noName));
end

% Disambiguate reserved names for table variables
if numNoName < numVars
    mdim = matlab.internal.tabular.private.metaDim;
    rlabels = mdim.labels;
    rnames = matlab.internal.tabular.private.varNamesDim.reservedNames;
    r = [rlabels rnames];
    varNames = matlab.lang.makeUniqueStrings(varNames,r,namelengthmax);
end
T = table.init(vars,m,{},numVars,varNames);
end