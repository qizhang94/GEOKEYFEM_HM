% For coupled flow and deformation problem, ndof + ndofp = 2 + 1 = 3
% sdof = (ndof + ndofp) * nnode

% A newer version with slightly different algorithms
% many_constant_matrices

function [K, all_sd_set_node, all_sd_B, all_sd_E] = pre_assemble_BigK(nod_adjele, area_nod, area_T3)

global ele_nods gcoord nnode nnel nel
global ndof ndofp

all_sd_set_node = cell(nnode, 1);
all_sd_B = cell(nnode, 1);
all_sd_E = cell(nnode, 1);

K11 = sparse(ndof*nnode, ndof*nnode);
K12 = sparse(ndof*nnode, ndofp*nnode);
K22 = sparse(ndofp*nnode, ndofp*nnode);

mat_Be = cell(nel, 1);
mat_Ee = cell(nel, 1);

for iel=1:nel       % loop for total number of element            
    for i=1:nnel    % loop for 3 nodes of (iel)-th element
        nod(i)=ele_nods(iel,i);  % extract nodes for (iel)-th element
        x(i)=gcoord(nod(i),1);   % extract x value of the node
        y(i)=gcoord(nod(i),2);   % extract y value of the node
    end
%     index=get_eledof(nod, nnel, ndofp);    % for flow
%     % Integration of N^T*N on one triangle element = constant matrix * area
%     % We can also add stablization term here (see follows)
%     % area_T3(iel)*[1/18, -1/36, -1/36; -1/36, 1/18, -1/36; -1/36, -1/36, 1/18];
%     K22(index, index) = K22(index, index) + area_T3(iel)*[1/6, 1/12, 1/12; 1/12, 1/6, 1/12; 1/12, 1/12, 1/6];
    
    %---------------------------------------
    % find the derivatives of shape functions
    %---------------------------------------	
    dNdx=(1/(2*area_T3(iel)))*[(y(2)-y(3)), (y(3)-y(1)), (y(1)-y(2))];  %derivatives w.r.t. x
    dNdy=(1/(2*area_T3(iel)))*[(x(3)-x(2)), (x(1)-x(3)), (x(2)-x(1))];  %derivatives w.r.t. y   
    Be = get_Bmat_2D(nnel,dNdx,dNdy);      % strain-displacement matrix of elements (for 2D CST3, it is a constant 3*6 matrix)                    
    mat_Be{iel}=Be;                      % save into matrix containing strain-displacement matrices of elements
    mat_Ee{iel} = [dNdx; dNdy];
end

for ino = 1:nnode    % Loop for the total number of node (= SD becuase of node-based)
    %----------------------------------------------------------------------
    % compute the required matrix of (ino)-th SD
    %----------------------------------------------------------------------
    ino_adjele = nod_adjele{ino};   % Matrix/row vector containing adjacent elements of each node
    ne = length(ino_adjele); % number of adjacent elements of node
    
    [nodB, B_ino, E_ino] = get_BEmat_NSFEM_T3_aver(nod_adjele{ino}, area_nod(ino), area_T3, mat_Be, mat_Ee);
    all_sd_set_node{ino} = nodB;
    all_sd_B{ino} = B_ino;
    all_sd_E{ino} = E_ino;
    
    index_u=get_eledof(nodB, length(nodB), ndof);    % extract system dofs for (ino)-th SD, for solid mechanics
    index_p=get_eledof(nodB, length(nodB), ndofp);    % extract system dofs for (ino)-th SD, for flow
    
    one_voigt = [1; 1; 0]; % Voigt notation of the second-order identity tensor
    K22_ino = E_ino'*E_ino*area_nod(ino);        % insert permeability/mobility tensor later

%     K11_ino=B_ino'*matmtx*B_ino*area_nod(ino);     % stiffness matrix of (ino)-th SD, for solid mechanics
%     K11(index_u, index_u) = K11(index_u, index_u) + K11_ino;
    
    K22(index_p, index_p) = K22(index_p, index_p) + K22_ino;
    
    % ino_adjele(ie)是当前focus on的element编号
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
        
        % 1*3 row vector [int_Np]
        [int_Np] = cal_int_N_T3(x_coord, y_coord, area_T3(ino_adjele(ie)), x_subSD, y_subSD);  % On one sub-SD, which is 1/3 of ino_adjele(ie)
        index_ele=get_eledof([n1, n2, n3], nnel, ndofp);
        K12(index_u, index_ele) = K12(index_u, index_ele) + B_ino'*one_voigt*int_Np;
    end
end

K = [K11, K12; transpose(K12), K22];     % Global tangential sparse matrix

end