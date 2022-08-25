function [STRESS, hsv, DDSDDE]=DP_UMAT(PROPS, STRESS0, DSTRAIN0, hsv0)
% Return mapping and CTO of the DP plasticity algorithm (2D plane strain)
% STRESS0: Old
% STRESS: New
% DDSDDE: Algorithmic consistent tangent operator
% Assume Ce (6*6) is uniform and constant

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR


E = PROPS(1); % Young's modulus
nu = PROPS(2); % Poisson's ratio
phi = PROPS(3);  %  phi  -  angle of friction (rad)
psi = PROPS(4);  %  psi  -  angle of dilation (rad)
coh = PROPS(5);  %  c  -  cohesion

mu = E/2/(1+nu);
K = E/3/(1-2*nu);

[Ce,~]=DlinElas(E,nu,6);

A = 3*sqrt(2)*coh*cos(phi)/sqrt(9+3*(sin(phi))^2);
B = 3*sqrt(2)*sin(phi)/sqrt(9+3*(sin(phi))^2);
b = 3*sqrt(2)*sin(psi)/sqrt(9+3*(sin(psi))^2);

stress_trial = STRESS0 + Ce*DSTRAIN0;   % Trial stress
stress_trial_tensor = [stress_trial(1), stress_trial(4), stress_trial(5); stress_trial(4), stress_trial(2), stress_trial(6); stress_trial(5), stress_trial(6), stress_trial(3)];
p_trial = 1/3*trace(stress_trial_tensor);
stress_dev_trial = stress_trial_tensor - p_trial*eye(3);
q_trial = sqrt(1.5)*norm(stress_dev_trial, 'fro');

if (sqrt(2/3)*q_trial + B*p_trial - A <= 0)
    % Not yield
    DDSDDE = Ce;
    STRESS = stress_trial;
    hsv = hsv0;
else
    % Yield, more complicated
    one_voigt = [1;1;1;0;0;0];
    I4 = diag([1,1,1,0.5,0.5,0.5]);
    nhat_tensor = stress_dev_trial/norm(stress_dev_trial, 'fro');
    nhat = [nhat_tensor(1,1); nhat_tensor(2,2); nhat_tensor(3,3); nhat_tensor(1,2); nhat_tensor(1,3); nhat_tensor(2,3)];
    delta_plastic = (sqrt(2/3)*q_trial + B*p_trial - A)/(2*mu + B*K*b);
    STRESS = stress_trial - delta_plastic*(K*b*one_voigt + 2*mu*nhat);
    hsv = hsv0 + sqrt(2/3)*delta_plastic; % deviatoric strain! Note that norm(nhat_tensor, 'fro') = 1
    DDSDDE = Ce - (K*b*one_voigt + 2*mu*nhat)*...
        transpose(K*B*one_voigt + 2*mu*nhat)/(2*mu + B*K*b)...
        - 4*mu^2*delta_plastic/norm(stress_dev_trial, 'fro')*(I4 - 1/3*(one_voigt*one_voigt') - nhat*nhat');
end

% See if DDSDDE is always invertible (it seems not!)
% coeff = 1/(sum(DDSDDE(1:3,1:3), 'all')/9);

% From 3D to 2D

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