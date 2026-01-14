function  [Js1,Js2] = MapsFacesCorners(DATA_REFMESH)

%%% Matrix that maps forces and moments defined over faces to forces and
%%% moments over corners
nmodesFACE = 6 ;
Js1 = cell(size(DATA_REFMESH)) ;
Js2 = cell(size(DATA_REFMESH)) ;

for itype = 1:length(DATA_REFMESH)
    As =DATA_REFMESH{itype}.Acondens_s1  ;
    J = (As') ;
    Js1{itype} = J(:,1:nmodesFACE) ;
    As =DATA_REFMESH{itype}.Acondens_s2  ;
    J =  (As') ;
    Js2{itype} = J(:,1:nmodesFACE) ;
end
