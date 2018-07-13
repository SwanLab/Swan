function [ DtJtil ] = cal_DtJtil_shear(smoothDtC,matCh,DtC,dim,ngaus)
% example: shear modulus maximization
% pp 747: S. amstutz et al.
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

n=2;
I = [1 3; 3 2];
DtC1 = zeros(npnod,1);
inv_C = inv(matCh);
DtJtil = [];

if (smoothDtC==1)
    DtC1 = zeros(npnod,1);
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    a=I(i,j);
                    b=I(k,l);
                    DtC1(:) = DtC(a,b,:);
                    DtJtil = DtJtil - 4*inv_C(3,a)*DtC1*inv_C(b,3)+2*inv_C(3,a)*DtC1*inv_C(b,3);
                end
            end
        end
    end
elseif (smoothDtC==0)
    DtC1 = zeros(ngaus,nelem);
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    a=I(i,j);
                    b=I(k,l);
                    DtC1(:,:) = DtC(a,b,:,:);
                    DtJtil = DtJtil - 4*inv_C(3,a)*DtC1*inv_C(b,3)+2*inv_C(3,a)*DtC1*inv_C(b,3);
                end
            end
        end
    end

end


















end

