function bc = bouncon(p,e,t,ldir,ldof,lval)

bc = [];

if ~isempty(ldir) & ~isempty(ldof) & ~isempty(lval)

    nlines = size(ldir,1);


    for j=1:nlines
        [ndir,~,~] = findnodes(p,e,t,ldir(j,:));
        
        if size(ndir,2)>1
            ndir = ndir';
        end
        
        
        ncc = size(ndir,1);

        node = [];
        dof = [];
        val = [];
        
        for i=1:ldof(j,1);

            node = [node;ndir];
            dof  = [dof; ldof(j,i+1)*ones(ncc,1)]; 
            val = [val;lval(j,i)*ones(ncc,1)];

        end

        bci = [node, dof, val];

        bc = [bc;bci];

    end
end

end