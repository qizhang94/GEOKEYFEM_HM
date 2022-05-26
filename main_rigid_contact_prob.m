% Rigid body contact with soft soil problem
restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel
global ndof ndofp bcdof bcval_incr g_const
format short;
addpath('.\SFEM_basic\');

%% Geometry
geoheight=1;  % meter
geowidth=1;
foot=0.15;

keypoint = [0, 0;
    geowidth, 0;
    geowidth, geoheight;
    0, geoheight;
    foot, geoheight;
    foot, 0;
    0, 0];
meshsize = [0.02, 0.02];
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

% MPa (E = 100 kPa)
K = 1/12; nu = 0.3; tau = 2; k_perm = 1e-14; mu_f = 1e-9; rho_buo = 0; % buoyant density of the mixture
data_const.tauG = tau/(3*K*(1 - 2*nu)/1/(1+nu));   % Stabilization multiplier, already divided by 2G
data_const.mobility = k_perm/mu_f;

cohesion = 0.005; phi = pi/12; psi = 0;
Props=ones(nnode, 1)*[3*K*(1 - 2*nu), nu, phi, psi, cohesion];

in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
[~, ~, cto_ela_t0] =...
        MohrCoulomb_UMAT(0, Props(1,:), in_situ_stress, zeros(3,1), 0);
data_const.Ce_t0 = cto_ela_t0;


old_solution = sparse((ndof + ndofp)*nnode, 1);
new_solution = sparse((ndof + ndofp)*nnode, 1);
residual_traction = sparse(ndof*nnode, 1); % not used

b_left=find(gcoord(:,1)==0);
b_bottom=find(gcoord(:,2)==0);
b_right=find(gcoord(:,1)==geowidth);
b_top=find(gcoord(:,2)==geoheight);
b_foot =intersect(find(gcoord(:,1)<=foot), b_top);


% Dirichlet BC
bcdof = [ndof*b_left'-1, ndof*b_right'-1, ndof*b_bottom', ndof*nnode + b_top'];
bcdof = unique(bcdof, 'stable');

% Neumann BC (do not need to do anything since we DO NOT prescribe traction)

num_steps = 10;
watch.dt = 8.64*ones(1,num_steps) ;   % time STEP INTERVAL, in second (= 0.0001 day), "pseudo" time
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
epsp = 1;

%% For contact model
rigid_left=[0, geoheight];
rigid_right=[foot, geoheight];
node_mas=[rigid_right; rigid_left]; % master node location

OMEGAN=1e8;
OMEGAT=1e7;
CFRI=0.3;
LTAN=1;
ct_prop=[OMEGAN,OMEGAT,CFRI,LTAN];

disptotal=-[0, 0.1];
disp_step=disptotal/num_steps;


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
stress = in_situ_stress*ones(1,nnode);  % Old
stress_new = stress;
SDV = zeros(1,nnode); % Solution-Dependent State Variables
SDV_new = SDV;

for step = 1:length(watch.dt)
    watch.now = watch.now + watch.dt(step);
    delta_t  = watch.dt(step);
    print_info = 'LOAD STEP = %d; TIME = %.2f:\n';
    fprintf(print_info, step, watch.now);

    for iter = 1:max_newton_iter

        [fCtotal_t,KCtotal_t]=fK_contact(new_solution, ct_prop, b_top, node_mas, disp_step);
        
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
        RHS(nnode*ndof+1:end) = RHS(nnode*ndof+1:end) - K22_correction*delta_t*new_solution(ndof*nnode+1:end);        

        % Contact mechanics
        RHS = RHS + fCtotal_t;
        K = K + KCtotal_t;

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

    node_mas = node_mas + [disp_step; disp_step];
    
end

clearvars K_example Mass_mat S_PPP int_Nu
clearvars ino step print_info delta_t error_temp sdof
clearvars iter
clearvars new_solution_temp

%% Plot and post-processing

seg = node_mas; % Final rigid body location
elemType = 'T3';

dispNodes = b_bottom;
dispNodes1=union(b_left, b_right);


u_xp = new_solution(1:2:nnode*ndof-1);
u_yp = new_solution(2:2:nnode*ndof);

node = gcoord;
element = ele_nods;

% Plot the FEM mesh 
figure('Color',[1 1 1])
hold on
plot_mesh(node,element,elemType,'k-');
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
title('Undeformed FE mesh')

% Plot numerical deformed configuration
fac=1;         
figure
clf
hold on
plot_mesh(node+fac*[u_xp u_yp],element,elemType,'k-');
title(' Numerical deformed mesh')
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);

% Plot deformation intensity with a colormap
figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_xp u_yp],element,elemType,u_xp);
colorbar
title('Deformation plot, U_X, m')

subplot(2,1,2);
plot_field(node+fac*[u_xp u_yp],element,elemType,u_yp);
colorbar
title('Deformation plot, U_Y, m')

% Plot stress with a colormap
% kPa
t1 = 1e3*stress(1,:);
t2 = 1e3*stress(2,:);
t3 = 1e3*stress(4,:);
t4 = 1e3*stress(3,:);

FontSize=16;
hh=figure('Name','Stress');
set(hh,'color','w');
set(hh,'Position',[400 100 1000 800]);
subplot(2,2,1);
snscontour(node, element,fac, u_xp,u_yp,t1','{\it\sigma}_{xx}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);

subplot(2,2,2);
snscontour(node, element,fac, u_xp,u_yp,t2','{\it\sigma}_{yy}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);

subplot(2,2,3);
snscontour(node, element,fac, u_xp,u_yp,t3','{\it\sigma}_{xy}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);

subplot(2,2,4);
snscontour(node, element,fac, u_xp,u_yp,t4','{\it\sigma}_{zz}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);
sgtitle('Stress tensor components, kPa');


% Plot equivlent plastic strain with a colormap
figure
t5=SDV;
snscontour(node, element,fac, u_xp,u_yp,t5','PEEQ');
title(' Equivalent plastic strain')
hold on
plot(seg(:,1),seg(:,2),'k-','linewidth',2);
