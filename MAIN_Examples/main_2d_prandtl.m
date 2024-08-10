%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified integration force is changed to 增量 (u and p)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Drucker-Prager model was adopted to approximate the Mohr-Coulomb model

restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel
global ndof ndofp bcdof bcval_incr g_const
format short;
addpath('./SFEM_basic/');

%% Geometry
geoheight=15;  % meter
geowidth=15;
foot=1;

keypoint = [0, 0;
    geowidth, 0;
    geowidth, geoheight;
    0, geoheight;
    foot, geoheight;
    foot, 0;
    0, 0];
meshsize = [3, 0.1];
[gcoord, ele_nods] = meshfooting('Non-uniform mesh', keypoint, meshsize, true); % This function will add more paths to the program

% TRR = ST3Element(ele_nods, gcoord);
% TRR.PlotElements;
% hold on

%% Material properties, boundary condition
nnode = size(gcoord, 1);
nnel = 3;
ndof = 2;
ndofp = 1;
nel = size(ele_nods, 1);
g_const = 0;       % set it to zero could eliminate the effect of gravity

% MPa (E = 1000 kPa)
K = 25/3; nu = 0.3; tau = 2; k_perm = 1e-14; mu_f = 1e-9; rho_buo = 0; % buoyant density of the mixture
data_const.tauG = tau/(3*K*(1 - 2*nu)/1/(1+nu));   % Stabilization multiplier, already divided by 2G
data_const.mobility = k_perm/mu_f;

cohesion = 0.001; phi = pi/36; psi = pi/36;
Props=ones(nnode, 1)*[3*K*(1 - 2*nu), nu, phi, psi, cohesion];

in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
[~, ~, cto_ela_t0] =...
        LinearElastic_UMAT(Props(1,:), in_situ_stress, zeros(3,1), 0);
data_const.Ce_t0 = cto_ela_t0;


old_solution = sparse((ndof + ndofp)*nnode, 1);
new_solution = sparse((ndof + ndofp)*nnode, 1);
residual_traction = sparse(ndof*nnode, 1); % not used

b_left=find(gcoord(:,1)==0);
b_bottom=find(gcoord(:,2)==0);
b_right=find(gcoord(:,1)==geowidth);
b_top=find(gcoord(:,2)==geoheight);
b_foot =intersect(find(gcoord(:,1)<=foot), b_top);
temp = [gcoord(b_foot, 1), b_foot];
temp = sortrows(temp);
b_foot = temp(:, 2);
clearvars temp


% Dirichlet BC
bcdof = [ndof*b_foot', ndof*b_left'-1, ndof*b_right'-1, ndof*b_bottom', ndof*nnode + b_top'];
bcdof = unique(bcdof, 'stable');

% Neumann BC
t_crit = 100; % seconds
F_surcharge = 0; % MPa
press_burden = 0; % MPa

% Not used here, so set it to 0, t_crit is useless
traction_f = @(x,t)(-(F_surcharge +...
    press_burden*(x<=foot).*(t > 0)*(t/t_crit)*(t<=t_crit) + press_burden*(x<=foot).*(t>t_crit)));  % "-" means compression

Total_U = -0.01;
num_steps = 100;
watch.dt = 40*ones(1,num_steps);   % time STEP INTERVAL (pseudo time here, 40 is arbitrary)
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
bcval_incr(1:num_steps, 1:length(b_foot)) = Total_U/num_steps;
epsp = 0.1;

%% For contact model



%% FEM Solving
% Find adjacent elements of each node
[nod_adjele] = get_nod_adjele;
% Compute the area of SD associated with node and element areas
[area_nod, area_T3] = cal_area_nod_T3(nod_adjele);

[K_example, all_sd_set_node, all_sd_B, all_sd_E] = pre_assemble_BigK(nod_adjele, area_nod, area_T3);
[Mass_mat, S_PPP, int_Nu] = pre_assemble_MassN(area_T3);
[K11_correction, K22_correction] = assemble_stab(nod_adjele, area_nod, area_T3,...
    all_sd_set_node, all_sd_B, all_sd_E,...
    data_const, epsp);

sdof = (ndof + ndofp)*nnode;
index_not_constrained = setdiff(1:sdof, bcdof);
max_newton_iter = 100; tol = 1e-3; tol2 = 1e-7;

cellstress = cell(length(watch.dt), 1); % Store stress tensor on every node for all time steps
cellUP = cell(length(watch.dt), 1); % Store the "new_solution" for all time steps
cellSDV = cell(length(watch.dt), 1);
stress = in_situ_stress*ones(1,nnode);  % Old
stress_new = stress;
SDV = zeros(1,nnode); % Solution-Dependent State Variables
SDV_new = SDV;

for step = 1:num_steps
    watch.now = watch.now + watch.dt(step);
    delta_t  = watch.dt(step);
    print_info = 'LOAD STEP = %d; TIME = %.2f:\n';
    fprintf(print_info, step, watch.now);
    residual_traction = assign_tractionBC2(residual_traction, b_top, traction_f, watch.now);

    for iter = 1:max_newton_iter

        
        [K, RHS, stress_new, SDV_new] = assemble_system(nod_adjele, area_nod, area_T3,...
        all_sd_set_node, all_sd_B, all_sd_E, ...
        new_solution, old_solution,...
        data_const, stress, SDV,...
        K_example, S_PPP, Mass_mat, int_Nu, residual_traction, delta_t, ...
        Props, false); % NS-FEM

        % NS-FEM to SNS-FEM
        K(1:nnode*ndof, 1:nnode*ndof) = K(1:nnode*ndof, 1:nnode*ndof) + K11_correction;
        K(nnode*ndof+1:end, nnode*ndof+1:end) = K(nnode*ndof+1:end, nnode*ndof+1:end) + K22_correction*delta_t;
        RHS(1:nnode*ndof) = RHS(1:nnode*ndof) - K11_correction*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode));
        RHS(nnode*ndof+1:end) = RHS(nnode*ndof+1:end) - K22_correction*delta_t*(new_solution(ndof*nnode+1:end) - old_solution(ndof*nnode+1:end));        


        if iter == 1
            % Adjust RHS according to bcval_incr, ONLY has an impact for non-zero Dirichlet BC
            RHS = RHS - K(:, bcdof)*bcval_incr(step,:)';
            new_solution(bcdof) = old_solution(bcdof) + bcval_incr(step,:)';
            error_temp = 1;
        end
        if (iter > 1) && (error_temp < tol)
            break;
        end

        new_solution_temp = new_solution; % Store the temporary value BEFORE Newton update

        % Solve and update
        new_solution(index_not_constrained) = new_solution(index_not_constrained)...
        + K(index_not_constrained, index_not_constrained)...
        \RHS(index_not_constrained);

        error_temp = norm(new_solution(1:nnode*ndof) - new_solution_temp(1:nnode*ndof))...
            /(norm(new_solution(1:nnode*ndof) - old_solution(1:nnode*ndof)) + eps);
        print_info = '\t NEWTON = %d; ERROR = %.4E; RES = %.4E\n';
        fprintf(print_info, iter, error_temp, norm(RHS(index_not_constrained)));
    end
    
    % Save state and store (converged) result
    stress = stress_new;
    old_solution = new_solution;
    SDV = SDV_new;
    cellstress{step} = stress;
    cellUP{step} = old_solution;
    cellSDV{step} = SDV;

    % Resistance force below the footing
    resi_force(step) = trapz(gcoord(b_foot, 1),stress(2, b_foot')');
    
end

clearvars K_example Mass_mat S_PPP int_Nu
clearvars ino step print_info delta_t error_temp sdof
clearvars iter
clearvars new_solution_temp

% Analytical is 6.49
figure(1);
set(gcf, 'Color', 'w');
plot(-Total_U/num_steps*(0:1:num_steps), [0,-resi_force]/cohesion/foot/1, 'b-d','Linewidth', 1, 'MarkerFaceColor','y'); grid on; % thickness in z direction: 1 unit
xlabel('Penetration depth (m)');
ylabel('Normalized Reaction force $q/c_{\rm u}$', 'Interpreter', 'LaTeX');

%% Plot and post-processing

elemType = 'T3';

dispNodes = b_bottom;
dispNodes1=union(b_left, b_right);


u_xp = new_solution(1:2:nnode*ndof-1);
u_yp = new_solution(2:2:nnode*ndof);
p_nodal = new_solution(nnode*ndof+1:end);

node = gcoord;
element = ele_nods;

% % Plot the FEM mesh 
% figure('Color',[1 1 1])
% hold on
% plot_mesh(node,element,elemType,'k-');
% plot(node(dispNodes,1),node(dispNodes,2),'ks');
% plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
% title('Undeformed FE mesh')

% % Plot numerical deformed configuration
fac=1;        
% figure
% clf
% hold on
% plot_mesh(node+fac*[u_xp u_yp],element,elemType,'k-');
% title(' Numerical deformed mesh')
% plot(node(dispNodes,1),node(dispNodes,2),'ks');
% plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
% hold on


% % Plot deformation intensity with a colormap
% figure
% clf
% subplot(2,1,1);
% plot_field(node+fac*[u_xp u_yp],element,elemType,u_xp);
% colorbar
% title('Deformation plot, U_X, m')

% subplot(2,1,2);
% plot_field(node+fac*[u_xp u_yp],element,elemType,u_yp);
% colorbar
% title('Deformation plot, U_Y, m')

% % Plot stress with a colormap
% % kPa
% t1 = 1e3*stress(1,:);
% t2 = 1e3*stress(2,:);
% t3 = 1e3*stress(4,:);
% t4 = 1e3*stress(3,:);
% 
% FontSize=16;
% hh=figure('Name','Stress');
% set(hh,'color','w');
% set(hh,'Position',[400 100 1000 800]);
% subplot(2,2,1);
% snscontour(node, element,fac, u_xp,u_yp,t1','{\it\sigma}_{xx}');
% set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
% hold on
% 
% 
% subplot(2,2,2);
% FontSize=16;
% snscontour(node, element,fac, u_xp,u_yp,t2','{\it\sigma}_{yy}');
% set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
% caxis([-7, 1]);
% hold on

% subplot(2,2,3);
% snscontour(node, element,fac, u_xp,u_yp,t3','{\it\sigma}_{xy}');
% set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
% hold on

% subplot(2,2,4);
% snscontour(node, element,fac, u_xp,u_yp,t4','{\it\sigma}_{zz}');
% set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
% hold on
% sgtitle('Stress tensor components, kPa');


% Plot equivalent plastic strain with a colormap
figure
t5=cellSDV{end};
snscontour(node, element,fac, u_xp,u_yp,t5','PEEQ');
title('Equivalent plastic strain');
caxis([0,0.02]);
hold on

% % Pore pressure
% figure
% snscontour(node, element,fac, u_xp,u_yp,p_nodal*1e3,'PP');
% title('Pore pressure, kPa')
% hold on

