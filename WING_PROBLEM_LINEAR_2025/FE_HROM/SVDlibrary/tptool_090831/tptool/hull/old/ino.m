function [ip,cs] = ino(U)
% INO - mátrix felbontás az INO feltétel szerint
%
% [IP,BASE]=INO(U)
% 
% Az U mátrix nxr-es, n>r és SN, azaz a sorok összege
% 1 kell legyen;
% 
% A kimenet az nxr-es IP és az rxr-es BASE mátrix, ahol
% a BASE SN, az IP pedig SN, NN(nemnegatív elemu) 
% és INO(minden oszlop legkisebb eleme 0
% 
% és amelyekre U=IP*BASE

%U=[.5 .1 .2 .2;.1 .1 .1 .7;.2 0 .4 .4; .2 .8 0 0;.3 .3 .2 .2];
[n r] = size(U);
UU = U;
U(:,r) = 0;

egy = min(U);  % a burkolószimplex lapjainak egyenlete/1
egy(r) = max(sum(U,2));
cs = ones(r,1)*egy(1:r-1); % a burkolószimplex csúcsai/1
for a = 1:r-1
    cs(a,a) = 2*egy(r) + egy(a) - sum(egy); % a burkolószimplexcsucsai/2
end
cs(:,r) = 1 - sum(cs,2); %az r. dimenzió
%for a=1:r
%    RNO(a,:)=sum(cs)-(r-1)*cs(a,:);% avégleges szimplex csúcsai
%end

%ip=UU*inv(RNO); %az interpolációs mátrix
ip = UU*inv(cs);  %az interpolációs mátrix
