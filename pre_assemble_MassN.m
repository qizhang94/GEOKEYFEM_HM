function [A, B, C] = pre_assemble_MassN(area_T3)

global ele_nods nnel nel nnode
global ndofp ndof

A = sparse(ndofp*nnode, ndofp*nnode);
B = sparse(ndofp*nnode, ndofp*nnode);   % Stabilization
C = sparse(2, ndof*nnode);   % for body force term

for iel=1:nel       % loop for total number of element            
    index=get_eledof(ele_nods(iel, :), nnel, ndofp);
    index2=get_eledof(ele_nods(iel, :), nnel, ndof);
    % Integration of N^T*N on one triangle element = constant matrix * area
    % area_T3(iel)*[1/6, 1/12, 1/12; 1/12, 1/6, 1/12; 1/12, 1/12, 1/6]
    % We can also add stablization term here (see follows)
    % area_T3(iel)*[1/18, -1/36, -1/36; -1/36, 1/18, -1/36; -1/36, -1/36, 1/18];
    A(index, index) = A(index, index) + area_T3(iel)*[1/6, 1/12, 1/12; 1/12, 1/6, 1/12; 1/12, 1/12, 1/6];
    B(index, index) = B(index, index) + area_T3(iel)*[1/18, -1/36, -1/36; -1/36, 1/18, -1/36; -1/36, -1/36, 1/18];
    C(:, index2) = C(:, index2) + area_T3(iel)*[1/3, 0, 1/3, 0, 1/3, 0; 0, 1/3, 0, 1/3, 0, 1/3];
end

end