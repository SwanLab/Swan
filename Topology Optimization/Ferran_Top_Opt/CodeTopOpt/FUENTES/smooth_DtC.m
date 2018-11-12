function [ nodDtC ] = smooth_DtC(DtC,nstre,npnod,ngaus,nelem,ndime,nnode,...
                    Msmooth,element,coordinatesn,coordinatesa,problembsc)

nodDtC = zeros(nstre,nstre,npnod);
gf = zeros(ngaus,nelem);
for i=1:nstre
    for j=1:nstre
        gf(:,:) = DtC(i,j,:,:);
        [gn] = smooth(nelem,npnod,ndime,nnode,gf,Msmooth,element,coordinatesn,coordinatesa,...
            problembsc);
        nodDtC(i,j,:)=gn;
    end
end
    

end

