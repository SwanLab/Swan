function [coord,connec,dim] = readFromGID(filename)
% *************************************************************************
%
%
% Read GID file (line-by-line)
%
%
% *************************************************************************
%   Inputs
%       - filename
%
%   Outputs
%       - coord: nodal coordinates
%       - connec: table of connectivities
%       - dim: dimension of the problem
%
% *************************************************************************

[path,name,ext] = fileparts(filename);
% Create a '.txt' file (fopen does not work w/ '.msh' files)
copyfile(filename, fullfile(path, [name '.txt']));
fid = fopen(filename,'r');
i = 1; % node index
j = 1; % elem index
while ~feof(fid)
    line = fgetl(fid); % Returns the next line
    data = strsplit(line,' ');
    first = data{1};
    
    % Avoid first-comment line
    comments = {'###','#'};
    if strcmp(first, comments) == 0
        
        % Find dimension
        if strcmp(first, 'MESH') == 1
            indice = find(strcmp(data,'dimension'));
            dim = str2double(data{indice+1});
            % Find nodal coordinates
        elseif strcmp(first, 'coordinates') == 1
            while strcmp(data{1}, 'end') == 0
                line = fgetl(fid);
                data = strsplit(line,' ');
                m = zeros(size(data,1),size(data,2));   % define dimensions
                m = str2double(data);                   % convert to double
                m = m(~isnan(m));                       % erase nan components
                % Save coordinates
                if isnumeric(m) == 1 && isempty(m) == 0
                    coord(i,:) = m;
                    i = i + 1;
                end
            end
            % Find connectivities between elements
        elseif strcmp(first, 'elements')
            while strcmp(data{1}, 'end') == 0
                line = fgetl(fid);
                data = strsplit(line,' ');
                m = zeros(size(data,1),size(data,2));   % define dimensions
                m = str2double(data);                   % convert to double
                m = m(~isnan(m));                       % erase nan components
                % Save connectivities
                if isnumeric(m) == 1 && isempty(m) == 0
                    connec(j,:) = m;
                    j = j + 1;
                end
            end
        end
    end
end
fclose(fid);
% Erase first column (NOT USED)
coord = coord(:,2:end);
connec = connec(:,2:end-1);
end

