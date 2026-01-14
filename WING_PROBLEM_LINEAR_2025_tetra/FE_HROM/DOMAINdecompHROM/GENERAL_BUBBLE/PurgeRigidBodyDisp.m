function SNAPcompl = PurgeRigidBodyDisp(DATAoffline,SNAPcompl,PhiRB,Mdom,MESH)
%--------------------------------------------------------------------------
% FUNCTION: PurgeRigidBodyDisp
%
% PURPOSE:
%   This function removes the rigid body (RB) components from a set of
%   displacement snapshots (`SNAPcompl`) by projecting them onto the
%   orthogonal complement of the subspace spanned by the rigid body modes
%   (`PhiRB`). The purging is applied according to the method specified in
%   the `DATAoffline` structure.
%
%   The rigid body projection operator is defined using a mass-orthogonal
%   projection (geometric mass matrix `Mdom`), which ensures the RB
%   components are filtered out in a physically consistent way.
%
% INPUTS:
%   - DATAoffline : Structure containing control parameters. It must include
%                   the field `PurgeRBody_Method`, which determines the
%                   strategy for purging rigid body modes. Valid options are:
%                     * 'INFINITESIMAL_ROTATIONS_AND_TRANSLATION' (default)
%                     * 'ONLY_TRANSLATIONS'
%                     * 'LARGE_ROTATIONS_AND_TRANSLATION' (not implemented)
%   - SNAPcompl   : Cell array of displacement snapshots to be purged of
%                   rigid body components. Each cell contains a matrix of
%                   size [nDOFs x nSnapshots].
%   - PhiRB       : Matrix whose columns span the space of rigid body modes
%                   (e.g., 3 for translations, 6 for 3D translations + rotations).
%
% OUTPUT:
%   - SNAPcompl   : The same cell array as input, but with each snapshot
%                   projected onto the deformation space (orthogonal to RB modes).
%
% DEPENDENCIES:
%   - SprojDEF_operator: Function that performs the projection onto the
%     orthogonal complement of the RB space, using the operator:
%
%         SprojDEF = I - PhiRB (PhiRBᵀ Mdom PhiRB)⁻¹ PhiRBᵀ Mdom
%
%     which filters out rigid body content from a given snapshot.
%
% REMARKS:
%   - The 'INFINITESIMAL_ROTATIONS_AND_TRANSLATION' option assumes
%     small-strain kinematics and is used by default. This may be revisited
%     for finite kinematics in future implementations.
%   - The 'LARGE_ROTATIONS_AND_TRANSLATION' mode is declared but not yet implemented.
%   - The projection tolerance is defined within `SprojDEF_operator`.
%
% AUTHOR:Joaquín A. Hernández,  15-May-2025, Comments by ChatGPT4


DATAoffline = DefaultField(DATAoffline,'PurgeRBody_Method','INFINITESIMAL_ROTATIONS_AND_TRANSLATION');

switch DATAoffline.PurgeRBody_Method
    case 'INFINITESIMAL_ROTATIONS_AND_TRANSLATION'
        disp('Purging of rigid body modes is still done assuming strains are small (23-Feb-2025)...')
        disp('This should be modified in future versions of the program')
        for iproj = 1:length(SNAPcompl)
            SNAPcompl{iproj} = SprojDEF_operator(PhiRB,Mdom,SNAPcompl{iproj}) ;
        end
    case 'ONLY_TRANSLATIONS'
        disp('This was tested out of curiosity, do not use it. ')
        if size(PhiRB,2) == 6
            PhiPURGE = PhiRB(:,1:3) ; 
        elseif  size(PhiRB,2) == 3
            PhiPURGE = PhiRB(:,1:2) ; 
        else
            error('Option not implemented')
        end
        
          for iproj = 1:length(SNAPcompl)
            SNAPcompl{iproj} = SprojDEF_operator(PhiPURGE,Mdom,SNAPcompl{iproj}) ;
          end
        
    case 'LARGE_ROTATIONS_AND_TRANSLATION'
      %   error('Option not implemented yet (15-May-2025)')
         
        SNAPcompl = PurgeRigidBodyDisp_LARGE(DATAoffline,SNAPcompl,Mdom,MESH)  ; 
        
        
        
end

