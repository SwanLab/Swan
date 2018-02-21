function compute_strain_fourier(dim,element,Bmat,d_u)

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;


% compute strain
strain = zeros(nstre,nelem);
for istre=1:nstre
    for inode=1:nnode
        for idime=1:nunkn
            ievab = nunkn*(inode-1)+idime;
            lnods = element.conectivities(:,inode);
%             strain1(1,istre,:)=strain1(1,istre,:)+Bmat(istre,ievab,:).*...
%                 (coorda(inode,idime,:)-coordn(inode,idime,:));
            
                 strain(istre,:)=strain(istre,:) +(squeeze(Bmat(istre,ievab,:)).*d_u(lnods))';
        end
    end
end



end