function  [V,F,nF,T,nT]=retry_remove_tetrahedrons(V,F,nF,T,nT)

% Make a list of tetrahedrons and outside faces to be removed, to 
% create a new outside boundary
remove_tetra=false(1,nT);
remove_face=false(1,nT);
for i=1:nF
    check=any(T(1:nT,:)==F(i,1),2)&any(T(1:nT,:)==F(i,2),2)&any(T(1:nT,:)==F(i,3),2);
    j=find(check,1,'first');
    if(~isempty(j)), 
        remove_tetra(j)=true;  remove_face(i)=true; 
    end
end
remove_tetra_list=unique(find(remove_tetra));
remove_face_list=find(remove_face);
F_rem=F(remove_face_list,:);


% Make a list of new faces 
nFnew=0; F_new=zeros(length(remove_tetra_list)*4,3);
for i=1:length(remove_tetra_list)
    j=remove_tetra_list(i);
    nFnew=nFnew+1; F_new(nFnew,:)=T(j,[1 2 3]);
    nFnew=nFnew+1; F_new(nFnew,:)=T(j,[1 4 2]); 
    nFnew=nFnew+1; F_new(nFnew,:)=T(j,[1 3 4]); 
    nFnew=nFnew+1; F_new(nFnew,:)=T(j,[3 2 4]); 
end            

% Only keep outside face, remove inside faces
remove_newface=false(1,nFnew);
for i=1:nFnew
    check=any(any(F_rem==F_new(i,1),2)&any(F_rem==F_new(i,2),2)&any(F_rem==F_new(i,3),2));
    if(check)
        remove_newface(i)=true;
    else
        check=any(F_new==F_new(i,1),2)&any(F_new==F_new(i,2),2)&any(F_new==F_new(i,3),2);
        j=find(check);
        if(length(j)==1)
        elseif(length(j)==2)
            remove_newface(j)=true;
        else
            error('not possible');
        end
    end
end
F_new(remove_newface,:)=[]; 

% Remove faces by replacing them
remove_face_list=sort(remove_face_list(:),1,'descend');
for i=1:length(remove_face_list), F(remove_face_list(i),:)=F(nF,:); nF=nF-1; end

% Remove tetrahedrons by replacing them
remove_tetra_list=sort(remove_tetra_list(:),1,'descend');
for i=1:length(remove_tetra_list), T(remove_tetra_list(i),:)=T(nT,:); nT=nT-1; end

% add the new faces
for i=1:size(F_new,1), nF=nF+1; F(nF,:)=F_new(i,:); end
