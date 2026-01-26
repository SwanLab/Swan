function [BasisINTFall,Indexes,TEXTP,BasisINTFcand,IndicesRB,Dcomp,nDOFsFACEall] ...
    = BasisINTFall_local(Vrb,FACES_GROUPS,BasisUdef,BasisRdef,M,DATAIN,...
    SinvVal_Udef,SinvVal_Rdef,TEXTP,MasterDOFS_perface,fI,BasisUrb,Mdom)

BasisINTFall = cell(size(Vrb)) ; % This cell array includes all possible interface modes
Indexes = cell(size(Vrb)) ;
for ifgroup=1:length(FACES_GROUPS)
    f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
    f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
    iface = FACES_GROUPS{ifgroup}(1) ;
    [BasisINTFall_ifacegroup,TEXTP,IndexesLOC ] = ...
        ReactionAndInterfaceLocalEner(BasisUdef,BasisRdef,f1,f2,...
        Vrb{iface},M{iface},DATAIN,SinvVal_Udef,SinvVal_Rdef,ifgroup,TEXTP,...
        MasterDOFS_perface{ifgroup})  ;
    BasisINTFall{FACES_GROUPS{ifgroup}(1)} = BasisINTFall_ifacegroup ;
    BasisINTFall{FACES_GROUPS{ifgroup}(2)} = BasisINTFall_ifacegroup ;
    Indexes{FACES_GROUPS{ifgroup}(1)} = IndexesLOC ;
    Indexes{FACES_GROUPS{ifgroup}(2)} = IndexesLOC ;
    
end


iacum = 0 ;
for iface =1:length(Indexes)
    fff = fieldnames(Indexes{iface}) ;
    for ifield = 1:length(fff)
        floc = fff{ifield} ;
        Indexes{iface}.(floc) =  Indexes{iface}.(floc) + iacum ;
    end
    iacum =   Indexes{iface}.DEFslave(end);
end


% ------------------------------------
% STEP 2. All candidates
% ----------------------
BasisINTFcand =  cell(size(BasisINTFall)) ;  % Sparse format
IndicesRB = [] ;
for iface = 1:length(BasisINTFall)
    IndicesRBloc =  Indexes{iface}.RB  ;
    IndicesRB = [IndicesRB,IndicesRBloc] ;
    BasisINTFcand{iface} = sparse([BasisINTFall{iface}]) ;
end

[nelems,nDOFsFACEall] = cellfun(@size,BasisINTFcand) ;

% IndicesDEF = 1:sum(nCAND) ;
% IndicesDEF(IndicesRB) = [] ;

% STEP 3. Form a diagonal basis matrix
% ------------------------------------
BasisINTFall = blkdiag(BasisINTFcand{:}) ;

% STEP 4. Matrix Dcomp
% -------------------------------------
Dcomp = Dmat_kinematicConstraint(fI,BasisUrb,BasisUdef,SinvVal_Udef,BasisINTFall,Mdom) ;