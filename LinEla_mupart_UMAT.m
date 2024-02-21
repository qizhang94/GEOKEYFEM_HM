function [STRESS, hsv, DDSDDE]=LinEla_mupart_UMAT(PROPS, STRESS0, DSTRAIN0, hsv0)

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR


E = PROPS(1); % Young's modulus
nu = PROPS(2); % Poisson's ratio

G = E/2/(1+nu);

Ce = G*diag([2,2,2,1,1,1]);
DDSDDE = Ce;

STRESS = STRESS0 + Ce*DSTRAIN0;

hsv = hsv0;

DDSDDE = [DDSDDE(1:2, 1:2), DDSDDE(1:2, 4); DDSDDE(4, 1:2), DDSDDE(4,4)];

end