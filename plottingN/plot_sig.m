function plot_sig(fac,u_x,u_y,elemType,stress)

% plots the color coded stress distribution in the finite element
% region along with color bar scale.
global node element
t1=zeros(size(element))';
t1(:,:)=stress(1,:,:);
t2=zeros(size(element))';
t2(:,:)=stress(2,:,:);
t3=zeros(size(element))';
t3(:,:)=stress(3,:,:);
t4=zeros(size(element))';
t4(:,:)=stress(4,:,:);

figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,t1');
colorbar
title('Stress plot, sigma_xx')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,t2');
colorbar
title('Stress plot, sigma_yy')

figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,t3');
colorbar
title('Stress plot, sigma_xy')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,t4');
colorbar
title('Stress plot, sigma_zz')
end % end of function



