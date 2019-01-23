function [DtJtil] = cal_DtJtil_bulk(smoothDtC,matCh,DtC,dim,ngaus)
% example: bulk modulus maximization
% pp 747: S. amstutz et al.

% smoothDtC: 0 o 1, 0: DtC is defined in gauss points
%                   1: DtC is defined on the nodes  
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

n=2;
I = [1 3; 3 2];
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
                    DtJtil = DtJtil - (inv_C(1,a)*DtC1*inv_C(b,1)+inv_C(1,a)*DtC1*inv_C(b,2)+...
                        inv_C(2,a)*DtC1*inv_C(b,1)+inv_C(2,a)*DtC1*inv_C(b,2));
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
                    DtJtil = DtJtil - (inv_C(1,a)*DtC1*inv_C(b,1)+inv_C(1,a)*DtC1*inv_C(b,2)+...
                        inv_C(2,a)*DtC1*inv_C(b,1)+inv_C(2,a)*DtC1*inv_C(b,2));
                end
            end
        end
    end

end

