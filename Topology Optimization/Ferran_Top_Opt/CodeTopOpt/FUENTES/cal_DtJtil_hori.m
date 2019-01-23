function [ DtJtil ] = cal_DtJtil_hori(smoothDtC,matCh,DtC,dim,ngaus)
% example: horizontal rigidity maximization
% pp 746: S. amstutz et al.
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

n=2;
I = [1 3; 3 2];
inv_C = inv(matCh);
DtJtil = [];

if (smoothDtC==1)
    % smoothing of DtC
    DtC1 = zeros(npnod,1);
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    a=I(i,j);
                    b=I(k,l);
                    DtC1(:,:) = DtC(a,b,:,:);
                    DtJtil = DtJtil - inv_C(1,a)*DtC1*inv_C(b,1);
                end
            end
        end
    end
else
    % DtC is defined at elemental level
    DtC1 = zeros(ngaus,nelem);
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    a=I(i,j);
                    b=I(k,l);
                    DtC1(:,:) = DtC(a,b,:,:);
                    DtJtil = DtJtil - inv_C(1,a)*DtC1*inv_C(b,1);
                end
            end
        end
    end
end

end

