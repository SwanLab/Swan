function [I,D] = rangeSearch(X,Y,r)
% Copyright (c) 2019, Matthieu Aussal, Ecole Polytechnique, CMAP      
% GNU General Public License v3.0. 
% Compute range search as matlab native function "rangesearch.m" from the
% Statistics_Toolbox

% Check toolbox licence
[lic,~] = license('checkout','Statistics_Toolbox');

% Native Mathworks
if lic
    [I,D] = rangesearch(X,Y,r);        
    
% Homemade 
else
    % Infos
    disp('No Statistics_Toolbox licence found, use homemade rangesearch')
    
    % Dimensions
    Nx = size(X,1);
    Ny = size(Y,1);
    
    % Positive values
    tmp = min([X;Y],[],1);
    X   = X - ones(Nx,1)*tmp;
    Y   = Y - ones(Ny,1)*tmp;
    
    % Boxes defined by 3 dimensional range indices
    ix         = floor(X./r) + 1;
    iy         = floor(Y./r) + 1;
    [idx,~,ix] = unique(ix,'rows');
    [jdx,~,iy] = unique(iy,'rows');
    
    % Split X inside boxes
    nbr    = accumarray(ix,1,[length(idx),1]);
    [~,ix] = sort(ix);
    Ix     = mat2cell(ix,nbr,1);
    
    % Split Y inside boxes
    nbr    = accumarray(iy,1,[length(jdx),1]);
    [~,iy] = sort(iy);
    Iy     = mat2cell(iy,nbr,1);
    
    % 3D translations indices
    [i,j,k] = ndgrid(-1:1,-1:1,-1:1);
    T       = [i(:) j(:) k(:)];
    
    % Boxes interaction
    Ixy = cell(27,1);
    for i = 1:27
        % Translations for X boxes
        tmp = idx + ones(size(idx,1),1) * T(i,:);
        
        % Box intersection between Y and X
        [Lia,Locb] = ismember(jdx,tmp,'rows');
        
        % Interactions among Y box
        Ixy{i} = [find(Lia) Locb(Lia)];
    end
    
    % Vectorial format and sorted rows among Y box
    Ixy = cell2mat(Ixy);
    Ixy = sortrows(Ixy);
    
    % Cell format
    nbr = accumarray(Ixy(:,1),1,[size(jdx,1) 1]);
    Ixy = mat2cell(Ixy(:,2),nbr,1);
    
    % Distance computation by tensor product :  |X-Y| = |X|^2 - 2*X.Y + |Y|^2
    XdotX = X(:,1).*X(:,1) + X(:,2).*X(:,2) + X(:,3).*X(:,3);
    YdotY = Y(:,1).*Y(:,1) + Y(:,2).*Y(:,2) + Y(:,3).*Y(:,3);
    A     = [XdotX      , -2*X , ones(Nx,1)];
    B     = [ones(Ny,1) ,   Y  , YdotY     ];
    
    % Output initialization
    I = cell(Ny,1);
    D = cell(Ny,1);
    
    % Loop among Y box
    for i = 1:size(Ixy,1)
        % X indices from X boxes
        ix = cell2mat(Ix(Ixy{i}));
        
        % Y indices from Y box
        iy = Iy{i};
        
        % Squared distance matrix
        Rxy2 = A(ix,:) * B(iy,:)';
        
        % Extract close data
        Rxy2 = abs((Rxy2<r^2) .* Rxy2);
        
        % Range search
        [idx,jdx,val] = find(Rxy2);
        
        % Column vector format
        if (size(idx,2) > 1)
            idx = idx';
            jdx = jdx';
            val = val';
        end
        
        % Final indices
        if ~isempty(idx)
            nbr   = accumarray(jdx,1,[length(iy) 1])';
            I(iy) = mat2cell(ix(idx)',1,nbr);
            D(iy) = mat2cell(sqrt(val)',1,nbr);
        end
    end
end
end


% % Graphical representation
% figure
% hold on
% for i = 1:length(Ix)
%     plot3(X(Ix{i},1),X(Ix{i},2),X(Ix{i},3),'*','Color',rand(1,3))
% %     pause
% end
% axis equal
% 
% % Graphical representation
% figure
% hold on
% for i = 1:length(Iy)
%     plot3(Y(Iy{i},1),Y(Iy{i},2),Y(Iy{i},3),'*','Color',rand(1,3))
% %     pause
% end
% axis equal