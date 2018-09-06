

function phishift = shift2n(direction,phi,conditions)
% SPACE SHIFT FUNCTION
% Using Neumann, or Neumann for the directional gradient boundary conditions, we just shift
% our matrix over one index in a certain direction
% in order to take derivatives. 
  % SHIFTS LEVEL SET FUNCTION WITH NEUMANN OR NEUMANN FOR THE GRADIENT CONDITIONS
  switch direction
     case 'w'           % SHIFT WEST
             [m,n] = size(phi) ;
             phishift(1:m,1:n-1) = phi(1:m,2:n) ;
             switch conditions
               case 'n'                
                    phishift(1:m,n) = phi(1:m,n) ;     % NEUMANN CONDITIONS
               case 'ng' 
                    phishift(1:m,n) = 2*phi(1:m,n)-phi(1:m,n-1) ; % NEUMANN FOR THE DIRECTIONAL GRADIENT CONDITIONS
             end                
     case 'e'           % SHIFT EAST
             [m,n] = size(phi) ;
             phishift(1:m,2:n) = phi(1:m,1:n-1) ;
             switch conditions
               case 'n' 
                    phishift(1:m,1) = phi(1:m,1) ;     % NEUMANN CONDITIONS
               case 'ng' 
                    phishift(1:m,1) = 2*phi(1:m,1)-phi(1:m,2) ;   % NEUMANN FOR THE DIRECTIONAL GRADIENT CONDITIONS
             end
     case 'n'           % SHIFT NORTH
             [m,n] = size(phi) ;
             phishift(1:m-1,1:n) = phi(2:m,1:n) ;
             switch conditions
               case 'n' 
                    phishift(m,1:n) = phi(m,1:n) ;      % NEUMANN CONDITIONS
               case 'ng' 
                    phishift(m,1:n) = 2*phi(m,1:n)-phi(m-1,1:n) ;  % NEUMANN FOR THE DIRECTIONAL GRADIENT CONDITIONS
             end  
     case 's'           % SHIFT SOUTH
             [m,n] = size(phi) ;
             phishift(2:m,1:n) = phi(1:m-1,1:n) ;
             switch conditions
               case 'n' 
                 phishift(1,1:n) = phi(1,1:n) ;      % NEUMANN CONDITIONS
               case 'ng'                 
                 phishift(1,1:n) = 2*phi(1,1:n)-phi(2,1:n) ;       % NEUMANN FOR THE DIRECTIONAL GRADIENT CONDITIONS
             end
      otherwise 
             error('SHIFT N,S,E, OR W?')
  end
end

