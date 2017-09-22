classdef Preprocess<handle
    % Only reads
    
    properties
    end
    
    methods(Static)
        function [coord,connec,dim] = readFromGiD(filename)
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
        
        function [fixnodes, forces] = getBC()
                        % Dirichlet
                        % Node - Dimension - Value
                        fixnodes = [
                            1 1 0
                            1 2 0
                            3 1 0
                            3 2 0
                            8 1 0
                            8 2 0
                            ];
            
                        % Neumann --> Fpunc (global)
                        forces = [
                            11 2 -1
                            ];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fixnodes = [
%                 1 1 0
%                 1 2 0
%                 3 1 0
%                 3 2 0
%                 8 1 0
%                 8 2 0
%                 18 1 0
%                 18 2 0
%                 29 1 0
%                 29 2 0
%                 45 1 0
%                 45 2 0
%                 62 1 0
%                 62 2 0
%                 83 1 0
%                 83 2 0
%                 109 1 0
%                 109 2 0
%                 136 1 0
%                 136 2 0
%                 168 1 0
%                 168 2 0
%                 201 1 0
%                 201 2 0
%                 234 1 0
%                 234 2 0
%                 279 1 0
%                 279 2 0
%                 321 1 0
%                 321 2 0
%                 370 1 0
%                 370 2 0
%                 418 1 0
%                 418 2 0
%                 470 1 0
%                 470 2 0
%                 525 1 0
%                 525 2 0
%                 582 1 0
%                 582 2 0
%                 649 1 0
%                 649 2 0
%                 713 1 0
%                 713 2 0
%                 783 1 0
%                 783 2 0
%                 852 1 0
%                 852 2 0
%                 923 1 0
%                 923 2 0
%                 1008 1 0
%                 1008 2 0
%                 1087 1 0
%                 1087 2 0
%                 1172 1 0
%                 1172 2 0
%                 1259 1 0
%                 1259 2 0
%                 1346 1 0
%                 1346 2 0
%                 1440 1 0
%                 1440 2 0
%                 1537 1 0
%                 1537 2 0
%                 1641 1 0
%                 1641 2 0
%                 1744 1 0
%                 1744 2 0
%                 1850 1 0
%                 1850 2 0
%                 1959 1 0
%                 1959 2 0
%                 2064 1 0
%                 2064 2 0
%                 2185 1 0
%                 2185 2 0
%                 2303 1 0
%                 2303 2 0
%                 2431 1 0
%                 2431 2 0
%                 2551 1 0
%                 2551 2 0
%                 2676 1 0
%                 2676 2 0
%                 2809 1 0
%                 2809 2 0
%                 2942 1 0
%                 2942 2 0
%                 3086 1 0
%                 3086 2 0
%                 3228 1 0
%                 3228 2 0
%                 3369 1 0
%                 3369 2 0
%                 3514 1 0
%                 3514 2 0
%                 3659 1 0
%                 3659 2 0
%                 3819 1 0
%                 3819 2 0
%                 3975 1 0
%                 3975 2 0
%                 ];
            
%             %% Force Prescribed
%             % Node                Dimension                Value
%             forces = [
%                 9915 1 0
%                 9915 2 -1
%                 ];
            
        end
        
    end
end



