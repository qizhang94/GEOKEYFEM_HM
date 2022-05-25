function [fCtotal,KCtotal]=fK_contact(new_solution, ct_prop, node_slv_set, node_mas, dUm)

global gcoord nnode ndof ndofp

OMEGAN=ct_prop(1);
OMEGAT=ct_prop(2);
CFRI=ct_prop(3);
LTAN=ct_prop(4);

dofs_all = nnode*(ndof+ndofp);
ux_t = new_solution(1:2:nnode*ndof-1);
uy_t = new_solution(2:2:nnode*ndof);
node = gcoord + [ux_t, uy_t]; % current configuration

KCtotal=sparse(dofs_all,dofs_all);
fCtotal=sparse(dofs_all,1);
for isl=1:length(node_slv_set) % set of slave nodes, a column vector
    ind=node_slv_set(isl);
    
    dof_ct_slave=get_eledof(ind, 1, ndof);
    posi_nd=node(ind,:)';
    
    ELXY(:,1)=posi_nd;
    ELXYP(:,1)=posi_nd;
    ELXY(:,2)=node_mas(1,:)'+dUm'; % master node
    ELXYP(:,2)=node_mas(1,:)';
    ELXY(:,3)=node_mas(2,:)'+dUm';
    ELXYP(:,3)=node_mas(2,:)';
    [FORCE, STIFF] = cntelm2d(OMEGAN, OMEGAT, CFRI, ELXY, ELXYP, LTAN);
    if not(isempty(FORCE))
        fCtotal(dof_ct_slave)=FORCE(1:2,1);
        KCtotal(dof_ct_slave,dof_ct_slave)=STIFF(1:2,1:2);
    end
    
end
end %element iteration end