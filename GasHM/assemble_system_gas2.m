% Excess pressure formulation, isotropic permeability
% NS-FEM
function [K, RHS, stress_new, SDV_new] = assemble_system_gas2(nod_adjele, area_nod, area_T3,...
    all_sd_set_node, all_sd_B, all_sd_E,...
    new_solution, old_solution,...
    material_data_const, stress, SDV,...
    K_example, S_PPP, ~, ~, residual_traction, delta_t,...
    Props, coupling)

global ndof ndofp nnode ele_nods gcoord nnel

K11 = sparse(ndof*nnode, ndof*nnode);
residual_1 = sparse(ndof*nnode, 1);
residual_2 = sparse(ndofp*nnode, 1);
K22_gascomp = sparse(ndofp*nnode, ndofp*nnode);

for ino = 1:nnode
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    index_p = get_eledof(nodB, length(nodB), ndofp); % Not used
    strain_ino_new = all_sd_B{ino}*new_solution(index_u);
    strain_ino_old = all_sd_B{ino}*old_solution(index_u);
    ino_adjele = nod_adjele{ino};   % Matrix/row vector containing adjacent elements of each node
    ne = length(ino_adjele); % number of adjacent elements of node    
    
    % Constitutive model
    [stress_new(:, ino), SDV_new(:, ino), cto] = DP_UMAT(Props(ino,:), stress(:, ino), strain_ino_new - strain_ino_old, SDV(:, ino));
    
    % Take care of the "sign" (+ or -)
    residual_1(index_u) = residual_1(index_u) - transpose(all_sd_B{ino})*stress_new([1,2,4], ino)*area_nod(ino);
    K11(index_u, index_u) = K11(index_u, index_u) + transpose(all_sd_B{ino})*cto*all_sd_B{ino}*area_nod(ino);

    %% Deal with gas compressibility induced new flux term (Nov 10, 2023)
    % all_sd_B: 3 rows; all_sd_E: 2 rows
    gradp_SD = all_sd_E{ino}*new_solution(nnode*ndof + index_p); % Gradient of gas pressure at one smoothing domain
    q_Darcy = -material_data_const.mobility*gradp_SD; % Darcy flux vector


    %% Deal with the gas compressibility term (Nov 9, 2023)
    % Because it need to be integrated on each sub-smoothing domain
    for ie = 1:ne    % loop for adjacent elements of node 'ino'
        n1 = ele_nods(ino_adjele(ie),1); % 1st node of (ie)-th adjacent element
        n2 = ele_nods(ino_adjele(ie),2); % 2nd node of (ie)-th adjacent element
        n3 = ele_nods(ino_adjele(ie),3); % 3rd node of (ie)-th adjacent element
        xx1 = gcoord(n1,1); yy1 = gcoord(n1,2); % coords of 1st node
        xx2 = gcoord(n2,1); yy2 = gcoord(n2,2); % coords of 2nd node
        xx3 = gcoord(n3,1); yy3 = gcoord(n3,2); % coords of 3rd node
        
        x_coord = [xx1; xx2; xx3];
        y_coord = [yy1; yy2; yy3];
        center_x = mean(x_coord);
        center_y = mean(y_coord);

        if ino == n1
            x_subSD = [xx1; (xx1 + xx2)/2; center_x; (xx1 + xx3)/2];
            y_subSD = [yy1; (yy1 + yy2)/2; center_y; (yy1 + yy3)/2];
        elseif ino == n2
            x_subSD = [xx2; (xx1 + xx2)/2; center_x; (xx2 + xx3)/2];
            y_subSD = [yy2; (yy1 + yy2)/2; center_y; (yy2 + yy3)/2];
        else
            x_subSD = [xx3; (xx1 + xx3)/2; center_x; (xx2 + xx3)/2];
            y_subSD = [yy3; (yy1 + yy3)/2; center_y; (yy2 + yy3)/2];
        end
        [x_subSD, y_subSD] = poly2ccw(x_subSD, y_subSD);   % Make sure to be counter-clockwise
        
        % 3*1 and 3*3 small matrices
        [subSD_RHS, subSD_K] = cal_int_gascompress(x_coord, y_coord, area_T3(ino_adjele(ie)), x_subSD, y_subSD);  % On one sub-SD, which is 1/3 of ino_adjele(ie)
        index_ele=get_eledof([n1, n2, n3], nnel, ndofp); % nnel = 3
        residual_2(index_ele) = residual_2(index_ele) + subSD_RHS;
        K22_gascomp(index_ele, index_ele) = K22_gascomp(index_ele, index_ele) + subSD_K;
    end
  
end
K12 = -K_example(1:ndof*nnode, ndof*nnode+1:end);   % K21 = transpose(K12) because of symmetry
residual_1 = residual_1 + (-K12)*new_solution(ndof*nnode+1:end);
% Add external force Eq. (33) of White(2008) paper to residual_1
residual_1 = residual_1 + residual_traction;

factor = material_data_const.tauG;
K22_uns = -K_example(ndof*nnode+1:end, ndof*nnode+1:end)*(material_data_const.mobility)*delta_t; % Note of time increment, K22 is (semi-)negative definite
K22 = K22_uns - factor*S_PPP;    % Pressure stabilization, factor = \tau/(2*G)
K22 = K22 - K22_gascomp; % Note the "minus" sign here!

% For no-flow BC, the Eq. (35) of White(2008) paper simply vanishes (multiply delta_t)
residual_2 = residual_2 + transpose(-K12)*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode)) + (-K22_uns)*new_solution(ndof*nnode+1:end);
residual_2 = residual_2 + factor*S_PPP*(new_solution(ndof*nnode+1:end) - old_solution(ndof*nnode+1:end));

if ~coupling
    K12(:)=0;
    residual_2(:)=0;
end

K = [K11, K12; transpose(K12), K22]; % not invertible WITHOUT prescribing BC
RHS = [residual_1; residual_2];


    function [subSD_RHS, subSD_K] = cal_int_gascompress(x_coord, y_coord, area, x_subSD, y_subSD)
            pgas_nodal = new_solution(nnode*ndof + [n1; n2; n3]); % 3*1 vector
            pgas_nodal_old = old_solution(nnode*ndof + [n1; n2; n3]); % 3*1 vector
            
            subSD_RHS = zeros(3,1);
            subSD_K = zeros(3,3);

            x_gauss = [-sqrt(3)/3; sqrt(3)/3; sqrt(3)/3; -sqrt(3)/3];
            y_gauss = [-sqrt(3)/3; -sqrt(3)/3; sqrt(3)/3; sqrt(3)/3];
            weight = ones(4,1);
            
            nnel_subSD = 4;  % Number of nodes per element ("element" is just a sub-smoothing domain)
            % We cannot name this "nnel_subSD" as "nnel", because it will
            % override the original global variable
            
            for i = 1:4   % here "4" is not nnel, but the number of Gauss integration points
                % One gauss integration point in the real domain
                [shapeq4, dNdrq4, dNdsq4]=get_shape_Q4(x_gauss(i),y_gauss(i)); % shapeq4 is a column vector
                detJ = det(cal_jacob_2D(nnel_subSD, dNdrq4, dNdsq4, x_subSD, y_subSD));   % determinant of the Jacobian mapping matrix for "subSD"
                
                % One gauss integration point in the real domain
                x = shapeq4' *  x_subSD;
                y = shapeq4' *  y_subSD;
            
                N1 = 1/(2*area)*(x_coord(2)*y_coord(3) - x_coord(3)*y_coord(2) + ...
                    (y_coord(2) - y_coord(3))*x + (x_coord(3) - x_coord(2))*y);
                N2 = 1/(2*area)*(x_coord(3)*y_coord(1) - x_coord(1)*y_coord(3) + ...
                    (y_coord(3) - y_coord(1))*x + (x_coord(1) - x_coord(3))*y);
                N3 = 1/(2*area)*(x_coord(1)*y_coord(2) - x_coord(2)*y_coord(1) + ...
                    (y_coord(1) - y_coord(2))*x + (x_coord(2) - x_coord(1))*y);
            

                pgas_Gauss = [N1, N2, N3]*pgas_nodal;
                pgas_Gauss_old = [N1, N2, N3]*pgas_nodal_old;
                K_f = pgas_Gauss; % Gas compressibility, same unit as "pressure"
                
                % Evaluate (do not forget the POROSITY here!)
                subSD_RHS = subSD_RHS + weight(i) * ([N1; N2; N3]*(pgas_Gauss - pgas_Gauss_old)*material_data_const.poro/K_f + ...
                    delta_t*[N1; N2; N3]*gradp_SD'*q_Darcy/K_f) * detJ;
                subSD_K = subSD_K + weight(i) * [N1; N2; N3]*pgas_Gauss_old*material_data_const.poro/K_f^2*[N1, N2, N3] * detJ;
            
            end
    end


end