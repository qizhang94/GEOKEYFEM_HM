function [Be_aug, Ee_aug] = expand_BE_SD(nod_T3, nodB, Be, Ee)
% For stabilization of nodal integration
% Expand the column of Be (add a lot of 0) such that Be has the same shape
% as B_ino

% nod_T3 = [n1, n2, n3]
% Be = mat_Be{ino_adjele(ie)}
% Ee = mat_Ee{ino_adjele(ie)}
% nodB  = nodes contributing (ino)-th smoothing domain

% integration area = 1/3*area_T3(ino_adjele(ie))

Be_aug = zeros(3, 2*length(nodB));
Ee_aug = zeros(2, length(nodB));

for i = 1:3
    index = find(nodB == nod_T3(i));
    Be_aug(:, 2*index-1:2*index) = Be(:, 2*i-1:2*i);
    Ee_aug(:, index) = Ee(:, i);
end

end