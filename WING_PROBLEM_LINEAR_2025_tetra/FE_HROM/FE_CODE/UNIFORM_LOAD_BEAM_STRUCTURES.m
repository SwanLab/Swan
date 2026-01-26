
%%%% DISTRIBUTED LOAD OVER SURFACES (per each slice)
% ---------------------------------
%dbstop('5')
INPUTS_LOC = DefaultField(INPUTS_LOC,'BEAMLOADS',[]) ;

INLOC = INPUTS_LOC.BEAMLOADS ;
INLOC = DefaultField(INLOC,'LOAD_UNIFORM',[]) ;

LOADglo = INLOC.LOAD_UNIFORM ;  % LOAD{iface}(idom) = [t_x,t_y,t_z], iface = 3, 4...  LAteral surfaces
%
ISLOCAL_DEF = cell(size(LOADglo)) ; % If is empty, or zero ---> Global coordinates
INPUTS_LOC.BEAMLOADS = DefaultField(INPUTS_LOC.BEAMLOADS,'ISLOCAL',ISLOCAL_DEF) ;
ISLOCAL = INPUTS_LOC.BEAMLOADS.ISLOCAL ;

CNb = cell(1,ndim) ; Tnod= cell(1,ndim) ;
nfaces = min(size(CONNECTb,2),length(LOADglo)) ;
for iface = 1:nfaces
    % Loop over faces of the SLICE 
    LOADface=  LOADglo{iface} ;
    if isempty(ISLOCAL{iface})
        ISLOCAL_face = zeros(size(LOADface,1),1) ;
    else
        ISLOCAL_face = ISLOCAL{iface} ;
    end
    % Number of surfaces of this type
    nslices = min(size(CONNECTb,1),size(LOADface,1) );
    LOADface = LOADface(1:nslices,:) ;
    
    for idom = 1:nslices
        % Applied force
        
        LOAD = LOADface(idom,:) ;
        
        % Corresponding connectivities
        CNbLOC = CONNECTb{idom,iface} ;
        
        if  (iface ==1 & idom ~= 1)   | (iface == 2 & idom ~=nslices)
            % Nothing done
        else
            
            
            if   ISLOCAL_face(idom) && any(LOAD)
                % Compute normals to each boundary element
                NORMALSv = NormalsBoundary(COOR,CNbLOC)  ;
                
                % Check that are all normals are oriented in the same
                % direction
                
                %    CNbLOC_localn = RenumberConnectivities(CNbLOC,1:size(NORMALSv,2)) ;
                % Force in the normal direction
                LOAD = LOAD(1)*ones(1,3) ; % To cheat the program
                % The input should be given in the first entry
                % And what about TLOC ? We have to arrange NORMALSv
                % somehow... But how ?
            end
            
            
            
            for idim = 1:length(CNb)
                tCOMP = LOAD(idim) ;
                if tCOMP ~= 0
                    CNb{idim} = [CNb{idim} ; CNbLOC] ;
                    if ISLOCAL_face(idom) == 0
                        TLOC = tCOMP*ones(size(CNbLOC)) ;
                    else
                        % Local coordinates  (force in the normal direction)
                        % -----------------
                        TLOC = tCOMP*repmat(NORMALSv(idim,:)',1,size(CNbLOC,2)) ;
                        %    TLOC =reshape(TLOC,size(CNbLOC,1),[]) ;
                    end
                    Tnod{idim} = [  Tnod{idim} ; TLOC]  ;
                end
            end
        end
    end
    
    
    
    
    
end


%%% POINT LOADS --------------
% ----------------------------



INLOC = INPUTS_LOC.BEAMLOADS ;

INLOC = DefaultField(INLOC,'GENERALIZED_FORCES_ENDS_BEAM',[]) ;

if  isempty(INLOC.GENERALIZED_FORCES_ENDS_BEAM)
    Fpnt =zeros(ndim*nnode,1) ; %
else
    FORCES = INLOC.GENERALIZED_FORCES_ENDS_BEAM ;
    Fpnt =  PointLoadsBeams(DOMAINVAR,COOR,CONNECTb,TypeElementB,FORCES,DATA)  ;
end





%
% if  ndim ==3
%
% Mt =  INPUTS_LOC.DISTLOAD.TORSION ;
% FACES = fieldnames(Mt) ;
% for iface = 1:length(FACES)
%     M = Mt.(FACES{iface}) ;
%     if  (M~=0)
%         [dimREF signREF] = Faces2Index(FACES{iface}) ;
%         % Identify nodes belonging to plane FACES{iface}
%         coorREF = min(signREF*COOR(BoundaryNodes,dimREF)) ; coorREF = signREF*coorREF(1) ;
%         rnodBASEloc = find(abs(COOR(BoundaryNodes,dimREF)-coorREF)<1e-10) ;
%         NODESb = BoundaryNodes(rnodBASEloc) ;
%         %% Identufy nodes of NODESb pertaining to line dimREF-dimCOMP
%         if dimREF == 1
%             dimCOMP = 3 ;  dimOTH = 2;
%         elseif dimREF == 2
%             dimCOMP = 3 ;   dimOTH =1;
%         else
%             dimCOMP = 1 ;  dimOTH =2;
%         end
%
%         % Linear distribution in the dimCOMP direction along axis dimOTH
%         % -----------------------------
%         yMIN =   min(COOR(NODESb,dimOTH)) ;
%         yMAX =   max(COOR(NODESb,dimOTH)) ;
%         y0 = 0.5*(yMIN+yMAX) ;
%         b = yMAX - yMIN ;
%         zMIN =   min(COOR(NODESb,dimCOMP)) ;
%         zMAX =   max(COOR(NODESb,dimCOMP)) ;
%         h = zMAX-zMIN ;
%         tauMAX = 6*M/h/b^2 ;
%      %   dbstop('58')
%         CONNECTloc = ElemBnd(CONNECTb,NODESb) ;
%         NODESbALL = CONNECTloc(:) ;
%         y = COOR(NODESbALL,dimOTH) ;
%         tau = tauMAX/(yMAX-y0)*(y-y0) ;
%         tau = reshape(tau',size(CONNECTloc,1),size(CONNECTloc,2)) ;
%
%         CNb{dimCOMP} = [CNb{dimCOMP} ; CONNECTloc] ;
%         Tnod{dimCOMP} = [  Tnod{dimCOMP} ; tau]  ;
%
%
%
%
%         %
%         %         %%% Upper and lower edges
%         %         signREFg = [-1 1] ;
%         %
%         %
%         %         for ilines = 1:2
%         %             signREF = signREFg(ilines) ;
%         %             coorREF = min(signREF*COOR(BoundaryNodes,dimCOMP)) ; coorREF = signREF*coorREF(1) ;
%         %             rnodBASEloc = find(abs(COOR(NODESb,dimCOMP)-coorREF)<1e-10) ;
%         %             NODESbLOC = NODESb(rnodBASEloc) ;
%         %             coorOTH = COOR(NODESbLOC,dimOTH) ;
%         %             cMED= 0.5*(max(coorOTH) + min(coorOTH)) ;
%         %             b =   (max(coorOTH) - min(coorOTH)) ;
%         %             % Half line
%         %             iP1 = find(coorOTH<=cMED) ;
%         %             NODESbL{1} = NODESb(iP1) ;
%         %             %   SIGNn(1) = 1;
%         %             % Other half line
%         %             iP2 = find(coorOTH>=cMED) ;
%         %             NODESbL{2} = NODESb(iP2) ;
%         %             %    SIGNn(2) = -1;
%         %
%         %             p = 4*M/b^2 ;
%         %
%         %             %   for ihalf = 1:2
%         %             DOFSloc = 3*(NODESbL{ilines}-1)  + dimCOMP ;
%         %             Fpnt(DOFSloc) = Fpnt(DOFSloc) + signREF*p*(b/2)/length(DOFSloc) ;
%         %             %  end
%         %
%         %
%         %         end
%
%
%     end
% end
%
% end
%




% INPUTS_LOC = DefaultField(INPUTS_LOC,'POINTLOADS',[]) ;
% INPUTS_LOC.POINTLOADS = DefaultField(INPUTS_LOC.POINTLOADS,'Fpoint',[]) ;
%
% if ~isempty(INPUTS_LOC.POINTLOADS)
% DOFS =  INPUTS_LOC.POINTLOADS.DOFS ;
% f =  INPUTS_LOC.POINTLOADS.Fpoint;
% Fpnt(DOFS) = f ;
% end
