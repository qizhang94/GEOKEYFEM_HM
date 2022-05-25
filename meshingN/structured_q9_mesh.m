function [node,element]=structured_q9_mesh(pt1,pt2,pt3,pt4,numx,numy)
%Forms the element and node matrices for Q9 element, with the element 
%and node matrices arranged in a nodal counterclockwise order
%(primary nodes-> secondary nodes then internal nodes)
nnx=numx+1;
nny=numy+1;
node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);

inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];

element=make_elem(node_pattern,numx,numy,inc_u,inc_v);

[element,node]=q4totq9(element,node,numx,numy);

end


