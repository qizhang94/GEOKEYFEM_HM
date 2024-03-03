function [STRESS, hsv, DDSDDE]=LinEla_mupart_UMAT(PROPS, STRESS0, DSTRAIN0, hsv0)

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR
%DSTRAIN0 = DSTRAIN0 - sum(DSTRAIN0(1:3))/3*[1;1;1;0;0;0];

E = PROPS(1); % Young's modulus
nu = PROPS(2); % Poisson's ratio

G = E/2/(1+nu);
one_voigt = [1,1,1,0,0,0]';
Ce = 2*G*(diag([1,1,1,0.5,0.5,0.5])-1/3*(one_voigt*one_voigt')); % change between 0 and 1
DDSDDE = Ce;

STRESS = STRESS0 + Ce*DSTRAIN0;

hsv = hsv0;

DDSDDE = [DDSDDE(1:2, 1:2), DDSDDE(1:2, 4); DDSDDE(4, 1:2), DDSDDE(4,4)];

end