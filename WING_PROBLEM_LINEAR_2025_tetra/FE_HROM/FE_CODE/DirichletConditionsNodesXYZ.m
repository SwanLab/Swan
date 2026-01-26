function [Gb,dR,DOFr,DOFm] = ...
    DirichletConditionsNodesXYZ(INPUTS_LOC,LABEL,NODES_ENTITIES,nnn,COOR,CONNECTb,TypeElementB,DATA)

% Loop over entities
G = {} ; uBAR  = {} ; DOFr = {} ; DOFm = {} ;
for ilines = 1:nnn    
    a = INPUTS_LOC.DISP.(LABEL){ilines} ;
    if size(COOR,2)==2
        TRANS = a(1:2);
        ROT = a(3) ;
    else
        TRANS = a(1:3);
        ROT = a(4:6) ;
    end
    NODES = NODES_ENTITIES{ilines} ;
    [sss,zzz] = cellfun(@size,a) ;
    Temp = cellfun(@isempty,TRANS) ;
    Remp = cellfun(@isempty,ROT) ;    
    
    if all(Temp==1) && all(Remp==1)
        % No constraint
        Gloc = [] ; uBARloc = [] ; DOFrLOC = [] ; DOFmLOC = [] ;
    elseif (all(Temp==0) && all(Remp==0)) ||  (all(Temp==0) && any(Remp==1))
        [Gloc,uBARloc,DOFrLOC,DOFmLOC] = ...
            PRESCRIBED_DISP_OVER_ENTITY(a,NODES,COOR,CONNECTb,TypeElementB,...
            DATA)  ;
    elseif any(Temp==1)
        % Rotations are not considered
        % ----------------------------
        DOFrLOC = [] ; DOFmLOC = [] ; Gloc = [] ; uBARloc = [];
        DOFloc = small2large(NODES,size(COOR,2)) ;
        for iii = 1:length(TRANS)
            if ~isempty(TRANS{iii})
                DOFrLOC = [DOFrLOC; [DOFloc(iii:size(COOR,2):end)]] ;
                uBARloc = [uBARloc; TRANS{iii}*ones(length(NODES),1)] ;
            end
        end
        
    end
    G{end+1} =  Gloc ;
    uBAR{end+1} = uBARloc ;
    DOFr{end+1} = DOFrLOC ;
    DOFm{end+1} = DOFmLOC;
    
    
    
end
DOFrMAT = cell2mat(DOFr') ;
[DOFrUNIQUE,IND]  = unique(DOFrMAT) ;

if length(DOFrUNIQUE) < length(DOFrMAT)
    
    DOFm = cell2mat(DOFm') ;
    
    if ~isempty(DOFm)
        error('Non-compatible option. Repeated DOFs when imposing BCs cannot exist with affine BCs')
    end
    
    uBAR = cell2mat(uBAR') ;
    
    
    DOFrMAT = DOFrUNIQUE;
    Gb = [] ;
    
    dR = uBAR(IND) ;
    DOFmMAT = [] ;
else
    DOFmMAT = cell2mat(DOFm') ;
    
    Gb = sparse(length(DOFrMAT),length(DOFmMAT)) ;
    dR = zeros(length(DOFrMAT),1) ;
    
    iini = 1 ;
    jini = 1;
    for iii = 1:length(DOFr)
        ifin = iini+ length(DOFr{iii})-1 ;
        jfin = jini+length(DOFm{iii})-1 ;
        if ~isempty(G{iii})
            Gb(iini:ifin,jini:jfin) = G{iii} ;
            
            
        end
        if ~isempty(uBAR{iii})
            dR(iini:ifin) = uBAR{iii}  ;
        end
        iini = ifin + 1;
        jini = jfin + 1;
        
    end   
    
end


DOFr = DOFrMAT ;
DOFm = DOFmMAT ;

