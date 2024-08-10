% Excess pressure formulation, isotropic permeability
% NS-FEM
% Be careful of the CONSTITUTIVE function and code variable: K22_uns!!!!!!!!
function [K, RHS, stress_new, SDV_new] = assemble_system(nod_adjele, area_nod, area_T3,...
    all_sd_set_node, all_sd_B, all_sd_E,...
    new_solution, old_solution,...
    material_data_const, stress, SDV,...
    K_example, S_PPP, Mass_mat, int_Nu, residual_traction, delta_t,...
    Props, coupling)

global ndof ndofp nnode

K11 = sparse(ndof*nnode, ndof*nnode);
residual_1 = sparse(ndof*nnode, 1);

for ino = 1:nnode
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    index_p = get_eledof(nodB, length(nodB), ndofp);
    strain_ino_new = all_sd_B{ino}*new_solution(index_u);
    strain_ino_old = all_sd_B{ino}*old_solution(index_u);
    
    % Constitutive model **************************************************
    [stress_new(:, ino), SDV_new(:, ino), cto] = DP_UMAT(Props(ino,:), stress(:, ino), strain_ino_new - strain_ino_old, SDV(:, ino));
    %[stress_new(:, ino), SDV_new(:, ino), cto] = LinEla_mupart_UMAT(Props(ino,:), stress(:, ino), strain_ino_new - strain_ino_old, SDV(:, ino));

    % Take care of the "sign" (+ or -)
    residual_1(index_u) = residual_1(index_u) - transpose(all_sd_B{ino})*stress_new([1,2,4], ino)*area_nod(ino);
    K11(index_u, index_u) = K11(index_u, index_u) + transpose(all_sd_B{ino})*cto*all_sd_B{ino}*area_nod(ino);
  
end
K12 = -K_example(1:ndof*nnode, ndof*nnode+1:end);   % K21 = transpose(K12) because of symmetry
residual_1 = residual_1 + (-K12)*new_solution(ndof*nnode+1:end);
% Add external force Eq. (33) of White(2008) paper to residual_1
residual_1 = residual_1 + residual_traction;

factor = material_data_const.tauG;
K22_uns = -K_example(ndof*nnode+1:end, ndof*nnode+1:end)*(material_data_const.mobility)*delta_t; % Note of time increment, K22 is (semi-)negative definite
% K22_uns = -Mass_mat*(material_data_const.mobility); % For nearly incompressibility u-p form
K22 = K22_uns - factor*S_PPP;    % Pressure stabilization, factor = \tau/(2*G)

% For no-flow BC, the Eq. (35) of White(2008) paper simply vanishes (multiply delta_t)
residual_2 = transpose(-K12)*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode)) + (-K22_uns)*new_solution(ndof*nnode+1:end);
residual_2 = residual_2 + factor*S_PPP*(new_solution(ndof*nnode+1:end) - old_solution(ndof*nnode+1:end));

if ~coupling
    K12(:)=0;
    residual_2(:)=0;
end

K = [K11, K12; transpose(K12), K22]; % not invertible WITHOUT prescribing BC
RHS = [residual_1; residual_2];


end