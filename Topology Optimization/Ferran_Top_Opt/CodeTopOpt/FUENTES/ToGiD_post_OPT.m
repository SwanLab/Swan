function ToGiD_post_OPT(dim,mstep,file_name,ngaus,element,tdisp,coordinatesn,gamma_gp,gamma_print,g_gp,phi,structural_values,problembsc)
                        
etype = element.type;
ptype = problembsc.problemtype;
ndime=dim.ndime; npnod=dim.npnod;nelem=dim.nelem; nnode=dim.nnode;
nunkn = dim.nunkn;
switch  etype
    case 'TRIANGLE'
        gtype = 'Triangle'; %gid type
    case 'QUAD'
        gtype = 'Quadrilateral';
    case 'HEXAHEDRA'
        gtype = 'Hexahedra';
end

% Leer tipo de post proceso (step,time)--------------
if isfield(problembsc,'kindpost')
    switch problembsc.kindpost
        case 'TIME'
            timestep = tp_time;
        case 'STEP'
            timestep = mstep;
        otherwise
            timestep = mstep;
            fprintf(1,'asignado postproceso "No. step".........\n')
    end
else
    timestep = mstep;   % por default no. de paso
end
% Escribe el fichero de resultados
res_file = strcat(file_name,'_',num2str(mstep),'.flavia.res');
fid = fopen(res_file,'w');

switch  etype
    case 'TRIANGLE'
        if nnode==3
            idxgp = [1 2 3]; job=2;
            gid_write_headerpost(fid,gtype,ngaus,job)
        elseif nnode==6
            idxgp = [];
        end
    case 'QUAD'
        if nnode==4
            idxgp = [1 2 3 4 ]; job =3;
            gid_write_headerpost(fid,gtype,ngaus,job)
        elseif nnode==8
            idxgp = [1 7 9 3 4 8 6 2 5]; job=1;
            gid_write_headerpost(fid,gtype,ngaus,job)
        end
end


switch problembsc.TYPE
    case 'MICRO'
        if (problembsc.ppinfo(1)==1)
            tstres = structural_values.tstres;
            vdisp = structural_values.vdisp;
            nbdat = size(tstres,1);
            ngaus = size(tstres,2);
            nstre = size(tstres,3);
            nelem = size(tstres,4);
            
            switch problembsc.phisical_type
                case {'ELASTIC'}
                    
                    for idat=1:nbdat
                        nameres = ['Stress' num2str(idat)];
                        stres = squeeze(tstres(idat,1:ngaus,1:nstre,1:nelem));
                        gid_write_gauss_tensorfield(fid,nameres,timestep,stres,idxgp,ngaus,nstre,nelem);
                    end
                    
                    if (problembsc.ppinfo(9)==1)
                        print_vdisp(fid,timestep,nunkn,npnod,nbdat,vdisp);
                    end
                    
                case {'THERMAL'}
                    for idat=1:nbdat
                        nameres = ['FLUX' num2str(idat)];
                        stres = shiftdim(tstres(idat,1:ngaus,1:nstre,1:nelem),1);
                        gid_write_gauss_vfield(fid,nameres,timestep,stres,idxgp)
                    end
                    
                    if (problembsc.ppinfo(9)==1)
                        print_disp(fid,timestep,nunkn,npnod,nbdata,vdisp);
                    end

            end
        end
        
    case 'MACRO'
        
        if (problembsc.ppinfo(1)==1)
            stres = structural_values.stres;
            strain = structural_values.strain;
            ngaus = size(stres,1);
            nstre = size(stres,2);
            nelem = size(stres,3);
            
            switch problembsc.phisical_type
                case {'ELASTIC'}
                    nameres = ['Stress'];
                    stres = squeeze(stres(1:ngaus,1:nstre,1:nelem));
                    gid_write_gauss_tensorfield(fid,nameres,timestep,stres,idxgp,ngaus,nstre,nelem);
                    
                    nameres = ['Strain'];
                    strain = squeeze(strain(1:ngaus,1:nstre,1:nelem));
                    gid_write_gauss_tensorfield(fid,nameres,timestep,strain,idxgp,ngaus,nstre,nelem);
                    
                    nameres = ['Energy'];
                    fobj = structural_values.fobj;
                    gid_write_gauss_sclfield(fid,nameres,timestep,fobj,idxgp);
             
                case {'THERMAL'}
                    nameres = ['FLUX'];
                    gid_write_gauss_vfield(fid,nameres,timestep,stres,idxgp)

            end
        end
end


if (problembsc.ppinfo(31)==1)
    nameres = 'Gamma FUNCTION NODAL ';
        gid_write_nodal_sclfield(fid,nameres,timestep,gamma_print);
    nameres = 'Gamma FUNCTION GAUSS ';
        gid_write_gauss_sclfield(fid,nameres,timestep,gamma_gp',idxgp);
    
end


if (problembsc.ppinfo(36)==1)
    nameres = 'Shear Modulus';
    mu = (structural_values.mu);    
%     if element.smoothing
%     gid_write_nodal_sclfield(fid,nameres,timestep,mu);
%     else
        gid_write_gauss_sclfield(fid,nameres,timestep,mu',idxgp);
%     end
    
    
    nameres = 'Bulk Modulus';
    kappa = (structural_values.kappa);
%     if element.smoothing
%     gid_write_nodal_sclfield(fid,nameres,timestep,kappa);
%     else
        gid_write_gauss_sclfield(fid,nameres,timestep,kappa',idxgp);
%     end
   
end

if (problembsc.ppinfo(9)==1)
    
    switch problembsc.phisical_type
        case {'ELASTIC'}
            nameres = 'DISPLACEMENT';
            gid_write_vfield(fid,nameres,timestep,tdisp);
        case {'THERMAL'}
            nameres = 'TEMPERATURE';
            gid_write_nodal_sclfield(fid,nameres,timestep,tdisp);
    end
    
end

if (problembsc.ppinfo(9)==1)
    nameres = 'F_ext';
    gid_write_vfield(fid,nameres,timestep,reshape(structural_values.fext,dim.ndime,[])');
end

% if (problembsc.ppinfo(31)==1)
%     nameres = 'g nodal';
%     gid_write_nodal_sclfield(fid,nameres,timestep,g_nodal);
% end

if (problembsc.ppinfo(33)==1)
    nameres = 'g_gp';
    gid_write_nodal_sclfield(fid,nameres,timestep,g_gp);
    
    nameres = 'phi';
    gid_write_nodal_sclfield(fid,nameres,timestep,phi)
    %gid_write_gauss_sclfield(fid,nameres,timestep,g_gp',idxgp);    
end


status = fclose(fid);

end
