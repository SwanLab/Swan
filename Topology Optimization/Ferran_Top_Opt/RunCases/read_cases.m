function problem2run = read_cases(fname)

if nargin == 0
    fname = 'cases2run.xlsx';
end

% Read data from excel
warning('Reading cases from file!');
fprintf('\nCHECK LIST\n');
fprintf('\t-Incremental optimization correct?\n');
fprintf('\t-Create optimization problem correct?\n')
fprintf('\t-Excel source data correct?\n');
fprintf('\t-MACRO/MICRO correct?\n');
fprintf('\t-Init design variable correct?\n');
fprintf('\n');
check = input('Continue(Y/N)? ','s');
if ~strcmpi(check,'y')
    error('User stop');
end

data = readtable(fname);
[~,nfields] = size(data); % number of fields and rows

% Complete missing data
for i = 1:nfields
    x = data.(i);
    
    % Check data type (numeric or cell)
    if isnumeric(x)
        non_valid = isnan(x);
        if ~non_valid
            data.(i) = x;
        else
            x = complete_data(x,non_valid);
            data.(i) = x;
        end
        
    elseif iscell(x)
        non_valid = cellfun('isempty', x);
        if ~non_valid
            data.(i) = x;
        else
            x = complete_data(x,non_valid);
            data.(i) = x;
        end
    else
        error('Data type not valid.');
    end
    
end

problem2run = data;
if isnumeric(problem2run.casenum) % convert to str casenum
    problem2run.casenum = num2cell(problem2run.casenum);
    fun = @(num) num2str(num);
    problem2run.casenum = cellfun(fun,problem2run.casenum,'UniformOutput',false);
end

fun = @(str,num) sprintf('%s\\',str); % sprintf('%s\\Case %s\\',str,num);
problem2run.save_directory = cellfun(fun,problem2run.save_directory,problem2run.casenum,'UniformOutput',false);

end

function x = complete_data (x,non_valid)

idx = find(~non_valid); % find indices of valid data
idx_non = find(non_valid); % find indices of invalid data

fill_with = zeros(length(idx_non),1);
for k = 1:length(idx_non)
    i = idx_non(k);
    aux = find(idx<i,1,'last'); % replace with first previous valid entry
    fill_with(k) = idx(aux);
end

x(idx_non) = x(fill_with);

end