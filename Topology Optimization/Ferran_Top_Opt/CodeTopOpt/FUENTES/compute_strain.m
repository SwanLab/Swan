<<<<<<< HEAD
function strain = compute_strain(dim,Bmat,element,coordinatesn,coordinatesa)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;


lnods = zeros(3,nelem);
coordn = zeros(3,2,nelem);
coorda = zeros(3,2,nelem);
for i=1:nnode
    lnods(i,:)= element.conectivities(:,i);
    for idime=1:ndime
        coordn(i,idime,:)= coordinatesn(lnods(i,:),idime);       
        coorda(i,idime,:)= coordinatesa(lnods(i,:),idime);
    end
end

% compute strain
strain = zeros(nstre,nelem);
for istre=1:nstre
    for inode=1:nnode
        for idime=1:nunkn
            ievab = nunkn*(inode-1)+idime;
            strain(istre,:)=strain(istre,:)+shiftdim(Bmat(istre,ievab,:).*...
                (coorda(inode,idime,:)-coordn(inode,idime,:)),1);
        end
    end
end

        switch element.material.subtype

            case 'PLANESTRES'
                epoiss = ones(1,nelem)*element.material.poiss;
                strain(nstre+1,:) = (-epoiss./(1-epoiss)).*(strain(1,:)+strain(2,:));

        end

=======
function strain = compute_strain(dim,Bmat,element,coordinatesn,coordinatesa)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;


dirichlet_data = zeros(3,nelem);
coordn = zeros(3,2,nelem);
coorda = zeros(3,2,nelem);
for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
    for idime=1:ndime
        coordn(i,idime,:)= coordinatesn(dirichlet_data(i,:),idime);       
        coorda(i,idime,:)= coordinatesa(dirichlet_data(i,:),idime);
    end
end

% compute strain
strain = zeros(nstre,nelem);
for istre=1:nstre
    for inode=1:nnode
        for idime=1:nunkn
            ievab = nunkn*(inode-1)+idime;
            strain(istre,:)=strain(istre,:)+shiftdim(Bmat(istre,ievab,:).*...
                (coorda(inode,idime,:)-coordn(inode,idime,:)),1);
        end
    end
end

        switch element.material.subtype

            case 'PLANESTRES'
                epoiss = ones(1,nelem)*element.material.poiss;
                strain(nstre+1,:) = (-epoiss./(1-epoiss)).*(strain(1,:)+strain(2,:));

        end

>>>>>>> refs/remotes/origin/master
end