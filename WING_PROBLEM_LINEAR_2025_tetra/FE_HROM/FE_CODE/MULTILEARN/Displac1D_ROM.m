function d = Displac1D_ROM(COORx,uBAR,BasisINTrb,BasisINT,MODES_RIGID,pINT,ndim)

dx  = zeros(size(COORx));
dy = dx ;
dz = dx ;
% Boundary conditions
%--------------------
% First node
dR_x = [] ;
dR_y = [] ;
dR_z = [] ;
r = [] ;
if ~isempty(uBAR{1,1})
    r =[1] ;
    dR_x = uBAR{1,1}(1) ;
    dR_y = uBAR{1,1}(2) ;
    if ndim == 3
    dR_z = uBAR{1,1}(3) ;
    
    end
end
if ~isempty(uBAR{end,3})
    r =[r length(dx)] ;
    dR_x = [dR_x; uBAR{end,3}(1)] ;
    dR_y = [dR_y; uBAR{end,3}(2)] ;
    if  ndim == 3
    dR_z = [dR_z; uBAR{end,3}(3)] ;
    end
end
l = 1:length(dx) ;
l(r) = [] ;
%%%%
dx(r) = dR_x ;
dy(r) = dR_y ;
if ndim == 3
dz(r) = dR_z ;
end


%% Basis Interface
NormBasisINTRB =   sqrt(sum(BasisINTrb.^2,1)) ;


%
nmodesINTF = size(BasisINT,2) ;
imode = 1;
III = find(MODES_RIGID == imode)  ;
if ~isempty(III)
    p_i =  pINT(III:nmodesINTF:end) ;
    dx(l) = p_i/NormBasisINTRB(imode) ;
end
%%%
imode = 2;
III = find(MODES_RIGID == imode)  ;
if ~isempty(III)
    p_i = pINT(III:nmodesINTF:end) ;
    dy(l) = p_i/NormBasisINTRB(imode) ;
end
%%%
if ndim == 3
imode = 3;
III = find(MODES_RIGID == imode)  ;
if ~isempty(III)
    p_i = pINT(III:nmodesINTF:end) ;
    dz(l) = p_i/NormBasisINTRB(imode) ;
end
end

if ndim == 3
d = [dx';dy';dz'] ;
else
   d = [dx';dy' ] ; 
end
d = d(:) ;