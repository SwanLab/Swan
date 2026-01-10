function [fextDOMglo,rRBglo,reactDOMrbGLO] = RigidBodyReactions(BasisUrb,nDOM,BasisRrb,FORCE_PROJECTS,INTENSITY_LOADS)


% % ---------------
nRB = size(BasisUrb,2) ; % number of rigid body modes
reactDOMrbGLO = cell(1,nDOM) ; % Cell containing the reactDOMrb variable
fextDOMglo = cell(1,nDOM) ;  % External forces applied on each subdomain
CovMixInv = inv(BasisUrb'*BasisRrb) ;
rRBglo = [] ;
for idom = 1:nDOM
    for iproj = 1:length(FORCE_PROJECTS)
        if isempty(fextDOMglo{idom})
            fextDOMglo{idom} = FORCE_PROJECTS{iproj}*INTENSITY_LOADS{iproj}(idom) ;
        else
            fextDOMglo{idom} =fextDOMglo{idom} +  FORCE_PROJECTS{iproj}*INTENSITY_LOADS{iproj}(idom) ;
        end
    end
    fextDOMrb = BasisUrb'*fextDOMglo{idom} ;  % Resultant (force x, force y and moment around reference point) of external forces
    rRB =- CovMixInv*fextDOMrb ;
    reactDOMrbGLO{idom} = BasisRrb*rRB ;
    rRBglo = [rRBglo;rRB] ;
end