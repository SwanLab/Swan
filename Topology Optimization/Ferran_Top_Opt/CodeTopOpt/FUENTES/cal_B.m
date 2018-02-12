function [ Bmat ] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype)

switch ptype
    case '2D'
        switch  etype
            case {'TRIANGLE','QUAD'}
                
                switch ftype
                    case {'ELASTIC'}
                        Bmat = zeros(nstre,nnode*nunkn,nelem);
                        for inode=1:nnode
                            j = nunkn*(inode-1)+1;
                            Bmat(1,j,:)=cartd(1,inode,:);
                            Bmat(2,j+1,:)=cartd(2,inode,:);
                            Bmat(3,j,:)=cartd(2,inode,:);
                            Bmat(3,j+1,:)=cartd(1,inode,:);
                        end
                        
                    case {'THERMAL'}
                        Bmat = zeros(nstre,nnode*nunkn,nelem);
                        for inode=1:nnode
                            j = nunkn*(inode-1)+1;
                            Bmat(1,j,:)=cartd(1,inode,:);
                            Bmat(2,j,:)=cartd(2,inode,:);
                        end
                end
                
            otherwise
                error('undifined element type')
        end
    case '3D'
       switch  etype
            case {'TETRAHEDRA','HEXAHEDRA'}
                nstre = 6;
                Bmat = zeros(nstre,nnode*nunkn,nelem);
                for inode=1:nnode
                    j = nunkn*(inode-1)+1;
                    % associated to normal strains 
                    Bmat(1,j,:)=cartd(1,inode,:);
                    Bmat(2,j+1,:)=cartd(2,inode,:);
                    Bmat(3,j+2,:)=cartd(3,inode,:);
                    % associated to shear strain, gamma12
                    Bmat(4,j,:)=cartd(2,inode,:);
                    Bmat(4,j+1,:)=cartd(1,inode,:);
                    % associated to shear strain, gamma13
                    Bmat(5,j,:)=cartd(3,inode,:);
                    Bmat(5,j+2,:)=cartd(1,inode,:);
                    % associated to shear strain, gamma23
                    Bmat(6,j+1,:)=cartd(3,inode,:);
                    Bmat(6,j+2,:)=cartd(2,inode,:);
                end
       end
end

end
