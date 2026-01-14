
%%%% DISTRIBUTED LOAD OVER SURFACES
% ---------------------------------
%dbstop('5')
t = INPUTS_LOC.DISTLOAD.t ;
FACES = fieldnames(t) ;
CNb = cell(1,ndim) ; Tnod= cell(1,ndim) ;
for iface = 1:length(FACES)
    tLOC = t.(FACES{iface}) ;
    if any(tLOC)
        [dimREF signREF] = Faces2Index(FACES{iface}) ;
        % Identify nodes belonging to plane FACES{iface}
        coorREF = min(signREF*COOR(BoundaryNodes,dimREF)) ; coorREF = signREF*coorREF(1) ;
        rnodBASEloc = find(abs(COOR(BoundaryNodes,dimREF)-coorREF)<TOL) ;
        NODESb = BoundaryNodes(rnodBASEloc) ;
        %%% Connectivity matrix
        CONNECTloc = ElemBnd(CONNECTb,NODESb) ;
        for idim = 1:length(CNb)
            tCOMP = tLOC(idim) ;
            if tCOMP ~= 0
                CNb{idim} = [CNb{idim} ; CONNECTloc] ;
                TLOC = tCOMP*ones(size(CONNECTloc)) ;
                Tnod{idim} = [  Tnod{idim} ; TLOC]  ;
            end
        end
    end
end

if  ndim ==3

Mt =  INPUTS_LOC.DISTLOAD.TORSION ;
FACES = fieldnames(Mt) ;
for iface = 1:length(FACES)
    M = Mt.(FACES{iface}) ;
    if  (M~=0)
        [dimREF signREF] = Faces2Index(FACES{iface}) ;
        % Identify nodes belonging to plane FACES{iface}
        coorREF = min(signREF*COOR(BoundaryNodes,dimREF)) ; coorREF = signREF*coorREF(1) ;
        rnodBASEloc = find(abs(COOR(BoundaryNodes,dimREF)-coorREF)<1e-10) ;
        NODESb = BoundaryNodes(rnodBASEloc) ;
        %% Identufy nodes of NODESb pertaining to line dimREF-dimCOMP
        if dimREF == 1
            dimCOMP = 3 ;  dimOTH = 2;
        elseif dimREF == 2
            dimCOMP = 3 ;   dimOTH =1;
        else
            dimCOMP = 1 ;  dimOTH =2;
        end
        
        % Linear distribution in the dimCOMP direction along axis dimOTH
        % -----------------------------
        yMIN =   min(COOR(NODESb,dimOTH)) ;
        yMAX =   max(COOR(NODESb,dimOTH)) ;
        y0 = 0.5*(yMIN+yMAX) ;
        b = yMAX - yMIN ;
        zMIN =   min(COOR(NODESb,dimCOMP)) ;
        zMAX =   max(COOR(NODESb,dimCOMP)) ;
        h = zMAX-zMIN ;
        tauMAX = 6*M/h/b^2 ;
     %   dbstop('58')
        CONNECTloc = ElemBnd(CONNECTb,NODESb) ;
        NODESbALL = CONNECTloc(:) ;
        y = COOR(NODESbALL,dimOTH) ;
        tau = tauMAX/(yMAX-y0)*(y-y0) ;
        tau = reshape(tau',size(CONNECTloc,1),size(CONNECTloc,2)) ;
        
        CNb{dimCOMP} = [CNb{dimCOMP} ; CONNECTloc] ;
        Tnod{dimCOMP} = [  Tnod{dimCOMP} ; tau]  ;
        
        
        
        
        %
        %         %%% Upper and lower edges
        %         signREFg = [-1 1] ;
        %
        %
        %         for ilines = 1:2
        %             signREF = signREFg(ilines) ;
        %             coorREF = min(signREF*COOR(BoundaryNodes,dimCOMP)) ; coorREF = signREF*coorREF(1) ;
        %             rnodBASEloc = find(abs(COOR(NODESb,dimCOMP)-coorREF)<1e-10) ;
        %             NODESbLOC = NODESb(rnodBASEloc) ;
        %             coorOTH = COOR(NODESbLOC,dimOTH) ;
        %             cMED= 0.5*(max(coorOTH) + min(coorOTH)) ;
        %             b =   (max(coorOTH) - min(coorOTH)) ;
        %             % Half line
        %             iP1 = find(coorOTH<=cMED) ;
        %             NODESbL{1} = NODESb(iP1) ;
        %             %   SIGNn(1) = 1;
        %             % Other half line
        %             iP2 = find(coorOTH>=cMED) ;
        %             NODESbL{2} = NODESb(iP2) ;
        %             %    SIGNn(2) = -1;
        %
        %             p = 4*M/b^2 ;
        %
        %             %   for ihalf = 1:2
        %             DOFSloc = 3*(NODESbL{ilines}-1)  + dimCOMP ;
        %             Fpnt(DOFSloc) = Fpnt(DOFSloc) + signREF*p*(b/2)/length(DOFSloc) ;
        %             %  end
        %
        %
        %         end
        
        
    end
end

end




%%% POINT LOADS --------------
% ----------------------------
% No point loads ....
Fpnt =zeros(ndim*nnode,1) ; %
INPUTS_LOC = DefaultField(INPUTS_LOC,'POINTLOADS',[]) ; 
INPUTS_LOC.POINTLOADS = DefaultField(INPUTS_LOC.POINTLOADS,'Fpoint',[]) ; 

if ~isempty(INPUTS_LOC.POINTLOADS) 
DOFS =  INPUTS_LOC.POINTLOADS.DOFS ;
f =  INPUTS_LOC.POINTLOADS.Fpoint;
Fpnt(DOFS) = f ; 
end
