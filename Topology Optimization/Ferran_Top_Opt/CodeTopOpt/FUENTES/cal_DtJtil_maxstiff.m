function [ DtJtil ] = cal_DtJtil_maxstiff(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem,fext)
% 

n=2;
I = [1 3; 3 2];
DtC1 = zeros(npnod,1);
strain0 = fext.macrostra;

hnorm = fext.h_C_0;

for a=1:3
    for b=1:3
        DtC1(:) = DtC(a,b,:);
        DtJtil = DtJtil - strain0(a)*DtC1(:)*strain0(b)/(strain0*strain0')/abs(hnorm);
    end
end


% sigma0 = matCh\strain0';
% 
% for a=1:3
%     for b=1:3
%         DtC1(:) = DtC(a,b,:);
%         DtJtil = DtJtil - sigma0(a)*DtC1(:)*sigma0(b)/(strain0*strain0');
%     end
% end
 

end

