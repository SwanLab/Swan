function meshData = readMshFile(filename)
    % This function reads a .msh file and extracts coordinates and elements
    % for the domains 'Sprout' and 'Gel'.
    
    % Initialize output structure
    meshData = struct();
    
    % Open the .msh file
    fid = fopen(filename, 'r');
    
    % Check if file opened successfully
    if fid == -1
        error('Could not open file.');
    end
    
    currentDomain = '';  % Track the current domain being processed
    
    % Read through the file line by line
    while ~feof(fid)
        tline = strtrim(fgetl(fid)); % Read a line and remove whitespace
        
        if startsWith(tline, 'MESH')
            % Extract domain name
            domainName = extractBetween(tline, '"', '"');
            currentDomain = domainName{1};
            meshData.(currentDomain) = struct('Coordinates', [], 'Elements', []);
        elseif contains(tline, 'Coordinates')
            coords = readCoordinates(fid);
            if ~isempty(coords)
                meshData.(currentDomain).Coordinates = coords;
            end
        elseif contains(tline, 'Elements')
            elements = readElements(fid);
            meshData.(currentDomain).Elements = elements;
        end
    end
    
    % Close the file
    fclose(fid);
end

function coords = readCoordinates(fid)
    % This function reads coordinates from the .msh file
    coords = [];
    while true
        tline = strtrim(fgetl(fid));
        if contains(tline, 'End Coordinates')
            break;
        end
        coords = [coords; str2num(tline)]; %#ok<*ST2NM>
    end
end

function elements = readElements(fid)
    % This function reads elements from the .msh file
    elements = [];
    while true
        tline = strtrim(fgetl(fid));
        if contains(tline, 'End Elements')
            break;
        end
        elements = [elements; str2num(tline)];
    end
end
