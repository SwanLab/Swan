function CONNECTb = RemoveREpeatedConnectivities(CONNECTb) 

%dbstop('4')
[CONNECTbSORT,iSORT ]= sort(CONNECTb,2) ; 
% Unique rows
[CONNECTbUNIQUE iB iC]= unique(CONNECTbSORT,'rows') ;  % Unique rows of CONNECTbSORT
iSORT = iSORT(iB,:) ; 

nnode = size(CONNECTbUNIQUE,2) ; nelem = size(CONNECTbUNIQUE,1) ; 

CONNECTbUNIQUE =  CONNECTbUNIQUE' ; 
CONNECTbUNIQUE = CONNECTbUNIQUE(:) ; 
CONNECTb = zeros(size(CONNECTbUNIQUE))  ; 
iSORT = iSORT' ; 
iSORT = iSORT(:) ; 

for inode = 1:nnode
    CONNECTb(inode:nnode:end) = CONNECTbUNIQUE(iSORT == inode) ; 
end

CONNECTb = reshape(CONNECTb',nnode,nelem) ; 
CONNECTb = CONNECTb' ; 