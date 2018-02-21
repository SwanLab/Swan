function [structural_values,vdisp,nbdata,post,fextpoint] = ...
    module_M(gamma,element,fixnodes,problembsc,coordinates,fext,dim)
% solve ku=f, compute topological derivative, and evaluate cost function, J

% basic dimensions
npnod=dim.npnod; nndof=dim.nndof; ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; neleq=dim.neleq; nunkn = dim.nunkn; nstre = dim.nstre;
% basic variables
[coordinatesn,coordinatesa] = init_coord(coordinates);
pressure = []; intvar_new = []; intvar_old = []; 
d_u = zeros(nndof,1);

% postprocess info
%post=problembsc.ppinfo;

% Compute free and fixed degree of freedom
[fix_df,free_df] = fix_free_degree_freedom(element.type,nndof,nunkn,fixnodes,problembsc);
% Impose displacements, update coordinates
switch problembsc.phisical_type
    case {'ELASTIC'}
        [coordinatesa] = actual_df('BC',free_df,coordinatesa,pressure,problembsc,element,fixnodes,d_u);
end
% compute stiffness magtrix
%[StifMat,post] = stiff_OPT(dim,element,problembsc,coordinatesn,coordinatesa,gamma);
[StifMat,post] = stiff_OPT_symmetric(dim,element,problembsc,coordinatesn,coordinatesa,gamma);

% point loads
[fextpoint] = ext_fnode(nndof,nunkn,fext.pointload,fext.fvol,problembsc.phisical_type);

switch problembsc.TYPE
    
    case 'MICRO'
        nbdata = nstre;
        vstrain=diag(ones(nbdata,1));
        
%         nbdata = nstre+1;
%         vstrain=[diag(ones(nbdata-1,1)); element.alpha'];
        
        
        % e11, e22 , 2e12
        [tstrain,tstres] = ini_tstrain_stres(dim,element,problembsc);
        structural_values.tstres = tstres;
        structural_values.tstrain = tstrain;
        structural_values.matCh = zeros(nstre,nstre);
        
    case 'MACRO'
        nbdata = 1;
        vstrain = zeros(nbdata,nstre);
        structural_values.costfunc = 0;
end

vdisp = zeros(nbdata,nndof);




for idata=1:nbdata
    
    if idata > nstre
        vstrain(idata,:) = (structural_values.matCh\element.alpha)';
    end
    strain_data = vstrain(idata,:); 
        
    [fextstra] = fext_strain(dim,element,problembsc,coordinatesn,coordinatesa,gamma,strain_data);
    fextlod = fextpoint + fextstra;
    
    % solve the system, ku=f
    [d_u,free1] = generalized_solver(d_u,StifMat,element,free_df,fix_df,dim,fixnodes,fextlod);
%     if (size(element.pnods,2)>0 || size(element.lglib,2)>0)
%         % periodic conditions
%         [d_u,free1] = solverp(d_u,StifMat,element.pnods,fextlod,fix_df,dim);
%     else
%         d_u = solver(d_u,free_df,StifMat,fextlod,fix_df,fixnodes);
%     end
    vdisp(idata,:)=d_u;
    % update coordinates
    
    switch problembsc.phisical_type
        case {'ELASTIC'}
            [coordinatesa] = actual_df('SOL',free_df,coordinatesn,pressure,problembsc,element,fixnodes,d_u);
    end
    
    
    % compute internal forces, force_int, and topological derivative, gfunc
    [force_int,structural_values] = cfint_OPT(dim,coordinatesn,coordinatesa,element,problembsc,gamma,idata,structural_values,vstrain,d_u,fextpoint);
    residual = force_int-fextlod;
    if (size(element.pnods,2)>0 || size(element.lglib,2)>0)
        % periodic conditions
        error = norm(residual(free1,1))/norm(fextlod);
    else
        error = norm(residual(free_df,1))/norm(fextlod);
    end
    
end


structural_values.d_u = d_u;
structural_values.fext = fextlod;
structural_values.vdisp = vdisp;
% switch problembsc.TYPE
%     case 'MACRO'
%         structural_values.costfunc = dot(fextlod,d_u);
% end


end



