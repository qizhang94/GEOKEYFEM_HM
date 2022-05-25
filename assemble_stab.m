function [K11_correction, K22_correction] = assemble_stab(nod_adjele, area_nod, area_T3,...
    all_sd_set_node, all_sd_B, all_sd_E,...
    material_data_const, epsp)

global ele_nods gcoord nnode nnel nel
global ndof ndofp g_const

K11_correction = sparse(ndof*nnode, ndof*nnode);   % For 积分稳定
K22_correction = sparse(ndofp*nnode, ndofp*nnode); % For 积分稳定

mat_Be = cell(nel, 1);
mat_Ee = cell(nel, 1);

for iel=1:nel       % loop for total number of element            
    for i=1:nnel    % loop for 3 nodes of (iel)-th element
        nod(i)=ele_nods(iel,i);  % extract nodes for (iel)-th element
        x(i)=gcoord(nod(i),1);   % extract x value of the node
        y(i)=gcoord(nod(i),2);   % extract y value of the node
    end
    dNdx=(1/(2*area_T3(iel)))*[(y(2)-y(3)), (y(3)-y(1)), (y(1)-y(2))];  %derivatives w.r.t. x
    dNdy=(1/(2*area_T3(iel)))*[(x(3)-x(2)), (x(1)-x(3)), (x(2)-x(1))];  %derivatives w.r.t. y   
    Be = get_Bmat_2D(nnel,dNdx,dNdy);      % strain-displacement matrix of elements (for 2D CST3, it is a constant 3*6 matrix)                    
    mat_Be{iel}=Be;                      % save into matrix containing strain-displacement matrices of elements
    mat_Ee{iel} = [dNdx; dNdy];
end

for ino = 1:nnode
    ino_adjele = nod_adjele{ino};   % Matrix/row vector containing adjacent elements of each node
    ne = length(ino_adjele); % number of adjacent elements of node
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    index_p = get_eledof(nodB, length(nodB), ndofp);
    for ie = 1:ne
        [Be_aug, Ee_aug] = expand_BE_SD(ele_nods(ino_adjele(ie), :), nodB, mat_Be{ino_adjele(ie)}, mat_Ee{ino_adjele(ie)});
        K11_correction(index_u, index_u) = K11_correction(index_u, index_u) + ...
            epsp*transpose(all_sd_B{ino} - Be_aug)*material_data_const.Ce_t0*(all_sd_B{ino} - Be_aug)*area_T3(ino_adjele(ie))/3;  % plus
        K22_correction(index_p, index_p) = K22_correction(index_p, index_p) - ...
            epsp*transpose(all_sd_E{ino} - Ee_aug)*(material_data_const.mobility)*(all_sd_E{ino} - Ee_aug)*area_T3(ino_adjele(ie))/3;   % minus
    end
end

end
