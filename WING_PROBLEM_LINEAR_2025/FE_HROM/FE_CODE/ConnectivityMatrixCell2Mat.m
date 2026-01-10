function CONNECTb = ConnectivityMatrixCell2Mat(CONNECTb)



if iscell(CONNECTb)
    % Sometimes this matrix is provided as a cell array (beam-like structures)
    % Here we convert it in a matrix
    CONNECTbMAT = [] ;
    % The boundary is formed by face 1 of domain 1, face 2 of domain N, and
    % the remaining faces of all domains
    
    if  length(size(CONNECTb)) == 2
        % CONNECTb is a cell array of N rows and M columns, N being the
        % number of slices in a beam-like structure, and M the number of
        % external surfaces. The 1st and 2nd of such
        %   surfaces corresponds  to the interface boundaries
        % ------------------------------------------------------------
        % JAHO, Nov-24th-2018. Something is wrong with this... Temporal
        % amendment
        if  size(CONNECTb,2) >1
            % Loop over external surfaces. The actual external boundary nodes are
            % those  belonging to face 1, slice 1, face 2, last slice, and
            % the remaining surfaces
            % ---------------------------
            for iface = 1:size(CONNECTb,2)
                if iface == 1
                    CONNECTbMAT = [CONNECTbMAT; (CONNECTb{1,iface})] ;
                elseif iface == 2
                    % Corrected 9-Dec-2018. This is only valid for beams.
                    % It may fe
                    CONNECTbMAT = [CONNECTbMAT; (CONNECTb{end,iface})] ;
                else
                    CONNECTbMAT = [CONNECTbMAT; cell2mat(CONNECTb(:,iface))] ;
                end
                
            end
        else
            CONNECTbMAT = cell2mat(CONNECTb) ;
            
        end
    else
        [ndomx,ndomy,nfaces ]= size(CONNECTb) ;
        for iface = 1:nfaces
            
            if iface == 1
                idomx = 1;
                idomy = 1:ndomy ;
            elseif iface == 3
                idomx = ndomx;
                idomy = 1:ndomy ;
            elseif iface == 2
                idomx = 1:ndomx;
                idomy = 1 ;
            elseif iface == 4
                idomx = 1:ndomx;
                idomy = ndomy ;
                
            else
                idomx = 1:ndomx ;
                idomy = 1:ndomy ;
            end
            
            
            CNLOC = CONNECTb(idomx, idomy  ,iface) ;
            CNLOC = cell2mat(CNLOC(:)) ;
            CONNECTbMAT = [CONNECTbMAT; CNLOC] ;
            
            
        end
        
        
    end
    
    
    
    CONNECTb = CONNECTbMAT ;
    
    
    
end

% Renumbering boundary connectivities
% (to ensure small bandwidth of both NstB and NstBw)

[~,ix]  = sort(CONNECTb(:,1)) ;
CONNECTb = CONNECTb(ix,:) ;