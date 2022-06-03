function [STRESS, hsv, DDSDDE]=LinearElastic_UMAT(PROPS, STRESS0, DSTRAIN0, hsv0)

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR


E = PROPS(1); % Young's modulus
nu = PROPS(2); % Poisson's ratio

[Ce,~]=DlinElas(E,nu,6);
DDSDDE = Ce;

STRESS = STRESS0 + Ce*DSTRAIN0;

hsv = hsv0;

DDSDDE = [DDSDDE(1:2, 1:2), DDSDDE(1:2, 4); DDSDDE(4, 1:2), DDSDDE(4,4)];

end


function [D,Dinv]=DlinElas(E,nu,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates linear elastic stiffness matrix D. If nsigma = 3, i.e.
% ! only three stress components plane stress is assumed
% !
% ! INPUT
% !     E   (real) Youngs module, module of elasticity
% !     nu    (real) Poisson's ratio
% !     nsigma  (int)  Number of stress components.
% !             nsigma = 4 in plane problems in general
% !             if nsigma = 3 plane stress is assumed
% !             nsigma = 6 in general stress states
% ! OUTPUT
% !     D   (real) Linear elastic stiffness matrix D
% !           Size is dependent on stress type
% !     Dinv  (real) Inverse of D, optional
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !------------------------------------------------------------------------

% implicit none
%
% integer(4), intent(in) :: nsigma
% real(8), intent(in) :: E, nu
%
% real(8), intent(out) :: D(nsigma,nsigma)
% real(8), intent(out), optional :: Dinv(nsigma,nsigma)
% !-------------------------
% real(8) c, g

D = zeros(nsigma,nsigma);
Dinv = zeros(nsigma,nsigma);
% !----- Plane stress with three stress components -----------------------
if (nsigma == 3)
    
    D(1,1) = 1.0 ; D(1,2) = nu;
    D(2,1) = nu ; D(2,2) = 1.0;
    D(3,3) = (1-nu)/2;
    D = E/(1-nu*nu) * D;
    
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0;
    Dinv(3,3) = 2.0*(1+nu);
    Dinv = Dinv/E;
    
    %         !----- Plane stress, strain or axisymmetry (four stress components) ---
elseif(nsigma == 4)
    D(1,1) = 1.0 - nu ; D(1,2) = nu         ; D(1,3) = nu;
    D(2,1) = nu         ; D(2,2) = 1.0 - nu ; D(2,3) = nu;
    D(3,1) = nu         ; D(3,2) = nu         ; D(3,3) = 1.0- nu;
    D(4,4) = (1-2.0*nu)/2;
    D = E/((1+nu)*(1-2.0*nu)) * D;
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu     ; Dinv(1,3) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0   ; Dinv(2,3) = -nu;
    Dinv(3,1) = -nu   ; Dinv(3,2) = -nu     ; Dinv(3,3) = 1.0;
    Dinv(4,4) = 2.0*(1+nu);
    Dinv = Dinv/E;
    %             !----- General stress state ---------------------------------------------
elseif(nsigma == 6)
    c = E/((1+nu)*(1-2.0*nu));
    g = E/(2.0*(1+nu));
    
    D(1,1) = (1-nu)*c ; D(1,2) = nu*c   ; D(1,3) = D(1,2);
    D(2,1) = D(1,2)   ; D(2,2) = D(1,1) ; D(2,3) = D(1,2);
    D(3,1) = D(1,2)   ; D(3,2) = D(1,2) ; D(3,3) = D(1,1);
    D(4,4) = g;
    D(5,5) = g;
    D(6,6) = g;
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu   ; Dinv(1,3) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0 ; Dinv(2,3) = -nu;
    Dinv(3,1) = -nu   ; Dinv(3,2) = -nu   ; Dinv(3,3) = 1.0;
    Dinv(4,4) = 2.0*(1+nu);
    Dinv(5,5) = 2.0*(1+nu);
    Dinv(6,6) = 2.0*(1+nu);
    Dinv = Dinv/E;
    
end


end %subroutine  DlinElas !-------------------------------------------------