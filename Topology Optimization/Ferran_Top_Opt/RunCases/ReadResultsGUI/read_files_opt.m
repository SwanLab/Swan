function [ output,message,all_iter_saved ] = read_files_opt( path )

%% Read files in path
files = dir(path);
files = struct2cell(files);
files = files(1,:)';

%% Read InfoCases
% Get file name
fname = 'InfoCases';
k = regexpi(files,fname);
index = ~cellfun('isempty',k);
if ~index
    message = 'Error: InfoCases file not found!';
    output = 0;
    all_iter_saved = 0;
    return;
end
fname = files{index};

% Open file and read
fid = fopen([path,'/',fname]);    

strpattern = 'Case %s -  %s %s %s %s - Vf %n Ptarget %n %*s %n %*s %n %*s %n %*s %n';
data = textscan(fid,strpattern);
fclose(fid);

% Export to table
varnames = {'NumCase','Case','Method','Algorithm','Kernel','Vtarget','Ptarget','lamV','penV','lamP','penP'};
output = table(data{:},'VariableNames',varnames);
    
%% Read RunTimes
% Get file name
fname = 'RunTimes';
k = regexpi(files,fname);
index = ~cellfun('isempty',k);
if ~index
    message = 'Error: RunTimes file not found!';
    output = 0;
    all_iter_saved = 0;
    return;
end
fname = files{index};

% Open file and read
N = length(output.NumCase);
output.RunTimes = cell(N,1);
fid = fopen([path,'/',fname]);
strpattern = 'Case %s - %n sec';
data = textscan(fid,strpattern);
fclose(fid);

% Export to table
output.RunTimes(1:length(data{2})) = num2cell(data{2});
output.RunTimes(cellfun('isempty',output.RunTimes)) = {'-'};

%% Read number of iterations
ncases = length(output.NumCase);
iter = zeros(ncases,1);
file = cell(ncases,1);
fig = cell(ncases,1);
all_iter_saved = cell(ncases,1);
for i = 1:ncases
    % Get last iteration file
    folder = sprintf('%s/Case %s - %s %s %s %s',path,output.NumCase{i},output.Case{i},output.Method{i},output.Algorithm{i},output.Kernel{i});
    if exist(folder,'dir')
        result_files = dir(folder);
        result_files = struct2cell(result_files);
        result_files = result_files(1,:)';
    else
        folder = sprintf('%s/Case %s -  %s %s %s %s',path,output.NumCase{i},output.Case{i},output.Method{i},output.Algorithm{i},output.Kernel{i}); % exception for compatibility
        if exist(folder,'dir')
            result_files = dir(folder);
            result_files = struct2cell(result_files);
            result_files = result_files(1,:)';
        else
            message = 'Error: Results directory not found!';
            output = 0;
            return;
        end
    end
    
    % Get number of iterations
    strpattern = [output.Case{i},'_%n.flavia.res'];
    fun = @(str) textscan(str,strpattern);
    aux = cellfun(fun,result_files);
    aux = cell2mat(aux);
    if isempty(aux)
        message = 'Error: Results file not found!';
        output = 0;
        return;
    else
        all_iter_saved{i} = unique(aux);
        iter(i) = max(aux); % total number of evaluations
    end
    
    % Get final results file
    file{i} = [folder,'/',output.Case{i},'_',num2str(iter(i)),'.flavia.res'];
    
    % Get results figure file
    fig{i} = [folder,'.fig'];
end
output.Iterations = iter;
output.ResultFile = file;
output.FigureFile = fig;
message = '';

end

