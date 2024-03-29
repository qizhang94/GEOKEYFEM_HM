%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified integration force is changed to "full" p and "incremental" u
% Could consider adsorption!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel
global ndof ndofp bcdof bcval_incr g_const
format short;
addpath('./SFEM_basic/');
addpath(genpath('./dengwirda-mesh2d-ceb68eb'));
addpath(genpath('./mesh2D_MLC'));
addpath('./meshingN');
addpath('./plottingN');

tic;

%% Geometry
node = [
    -1., -1.; +1., -1.
    +1., +1.; -1., +1.
    ];
edge = [
    1 ,  2 ;  2 ,  3
    3 ,  4 ;  4 ,  1
    ];

adel = 2.*pi / +32;
amin = 0.*pi;
amax = 2.*pi - adel;

xcir = +.1 * ...
    cos(amin:adel:amax)';
ycir = +.1 * ...
    sin(amin:adel:amax)';
ncir = [xcir,ycir];
numc = size(ncir,1);

ecir(:,1) = ...
    [(1:numc-1)'; numc];
ecir(:,2) = ...
    [(2:numc-0)'; +1  ];

ecir = ecir+size(node,1);
edge = [edge; ecir];
node = [node; ncir];

%---------------------------------------------- do mesh-gen.
hfun = @hfun8;

[vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun);

%---------------------------------------------- do mesh-opt.
[vert,~, ...
    tria,~] = smooth2(vert,etri,tria,tnum);

% figure;
% set(gcf,'Color','w');
% patch('faces',tria(:,1:3),'vertices',vert, ...
%     'facecolor','w', ...
%     'edgecolor',[.2,.2,.2]);
% hold on; axis image off;
% patch('faces',edge(:,1:2),'vertices',node, ...
%     'facecolor','w', ...
%     'edgecolor',[.1,.1,.1], ...
%     'linewidth',1.5);
% title(['MESH-OPT.: KIND=DELFRONT, $\rm \|TRIA\|$=', ...
%     num2str(size(tria,1))],'Interpreter','latex','fontsize',13);
% exportgraphics(gcf,'./gas_data/FE_mesh.pdf','ContentType','vector');

gcoord = vert; ele_nods = tria;

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

% MPa (following K value is in fact the generalized bulk modulus)
K = 2.8662e+03; nu = 0.15; tau = 0; k_perm = 2e-19; mu_f = 2e-11; rho_buo = 0; poro = 0.06; % buoyant density of the mixture
data_const.tauG = tau/(3*K*(1 - 2*nu)/1/(1+nu));   % PPP Stabilization multiplier, already divided by 2G
data_const.mobility = k_perm/mu_f;
data_const.poro = poro; % Porosity appears in the gas compressibility term
data_const.rhorock = 0; % kg/m^3
data_const.V_L = 0.015; % m^3/kg
data_const.P_L = 4; % MPa
data_const.T = 353.15; % Reservoir temperature [K]

% cohesion = 100; phi = pi/3; psi = pi/3; % make it elastic
Props = ones(nnode, 1)*[1.07, 0.01, 0, pi/6, poro]; % pi/6 is the bedding orientation

in_situ_stress = [-5; -10; -5; 0; 0; 0]; % initial effective stress field, in general, the first three components could be negative
Pc0 = -20; % Initial preconsolidation pressure, MPa

solid_data = Emptymodel(in_situ_stress, Pc0);
[~, ~, ~, cto_ela_t0] = AMCC_UMAT(Props(1,:), solid_data.stress_old, zeros(3,1), solid_data.hsv_old);
data_const.Ce_t0 = cto_ela_t0;

b_left=find(gcoord(:,1)== -1);
b_bottom=find(gcoord(:,2)== -1);
b_right=find(gcoord(:,1)== 1);
b_top=find(gcoord(:,2)== 1); % Prescribed gas pressure (through the penalty method)
% b_foot =intersect(find(gcoord(:,1)<=foot), b_top);

b_well = (5:36)'; % Inner boundary

old_solution = sparse((ndof + ndofp)*nnode, 1);
old_solution(nnode*ndof+1:end) = 30; % Initial uniform gas pressure of 2 [MPa]
old_solution(nnode*ndof+b_well) = 2; % Dirichelet boundary condition of p_gas
new_solution = old_solution;
traction_temp = sparse(ndof*nnode, 1);


% Dirichlet BC
bcdof = [ndof*b_well'-1, ndof*b_well', nnode*ndof+b_well'];
bcdof = unique(bcdof, 'stable');

% Neumann BC
F_x = 35; F_y = 40; 
traction_xL = @(x,t)(F_x); traction_xR = @(x,t)(-F_x); 
traction_yB = @(x,t)(F_y); traction_yT = @(x,t)(-F_y);


watch.dt = [10*ones(1,10), 10*(1.5.^(1:1:30)), 30*24*3600*ones(1,15)];   % time STEP INTERVAL
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
epsp = 1;
watch.max = 500*24*3600; % Maximum simulation time 500 days

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
SDV = zeros(2,nnode); % Solution-Dependent State Variables
SDV(1,:) = Pc0;
SDV_new = SDV;

for step = 1:length(watch.dt)
    watch.now = watch.now + watch.dt(step);
    delta_t  = watch.dt(step);
    if watch.now > watch.max
        delta_t = delta_t - (watch.now - watch.max);
        watch.now = watch.max;
    end     
    print_info = 'LOAD STEP = %d; TIME = %.2f:\n';
    fprintf(print_info, step, watch.now);
    residual_traction = assign_tractionBC2(traction_temp, b_right, traction_xR, watch.now, 'vertical') + ...
        assign_tractionBC2(traction_temp, b_left, traction_xL, watch.now, 'vertical') + ...
        assign_tractionBC2(traction_temp, b_top, traction_yT, watch.now) + ...
        assign_tractionBC2(traction_temp, b_bottom, traction_yB, watch.now);

    for iter = 1:max_newton_iter
        
        [K, RHS, stress_new, SDV_new] = assemble_system_gas(nod_adjele, area_nod, area_T3,...
        all_sd_set_node, all_sd_B, all_sd_E, ...
        new_solution, old_solution,...
        data_const, stress, SDV,...
        K_example, S_PPP, Mass_mat, int_Nu, residual_traction, delta_t, ...
        Props, true); % NS-FEM

        % NS-FEM to SNS-FEM
        K(1:nnode*ndof, 1:nnode*ndof) = K(1:nnode*ndof, 1:nnode*ndof) + K11_correction;
        K(nnode*ndof+1:end, nnode*ndof+1:end) = K(nnode*ndof+1:end, nnode*ndof+1:end) + K22_correction*delta_t;
        RHS(1:nnode*ndof) = RHS(1:nnode*ndof) - K11_correction*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode));
        RHS(nnode*ndof+1:end) = RHS(nnode*ndof+1:end) - K22_correction*delta_t*(new_solution(ndof*nnode+1:end) - 0*old_solution(ndof*nnode+1:end));        


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
    
end

clearvars K_example Mass_mat S_PPP int_Nu
clearvars ino step print_info delta_t error_temp sdof
clearvars iter
clearvars new_solution_temp

%% Plot and post-processing

toc;

elemType = 'T3';

dispNodes = b_bottom;
dispNodes1=union(b_left, b_right);

plotstep = 40; % 66.5810 days

u_xp = cellUP{plotstep}(1:2:nnode*ndof-1);
u_yp = cellUP{plotstep}(2:2:nnode*ndof);
p_nodal = cellUP{plotstep}(nnode*ndof+1:end);

node = gcoord;
element = ele_nods;

% Plot the FEM mesh 
figure('Color',[1 1 1])
hold on
plot_mesh(node,element,elemType,'k-');
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
title('Undeformed FE mesh');

% Plot numerical deformed configuration
fac=1;         
figure
clf
hold on
plot_mesh(node+fac*[u_xp u_yp],element,elemType,'k-');
title(' Numerical deformed mesh');
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
hold on


% Plot deformation intensity with a colormap
figure
clf
set(gcf,'Color','w');
subplot(2,1,1);
plot_field(node+fac*[u_xp u_yp],element,elemType,u_xp);
c = colorbar; c.FontSize = 13; set(c,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',15);
title('Deformation plot, $U_X$, m','Interpreter','latex','fontsize',15);

subplot(2,1,2);
plot_field(node+fac*[u_xp u_yp],element,elemType,u_yp);
c = colorbar; c.FontSize = 13; set(c,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',15);
title('Deformation plot, $U_Y$, m','Interpreter','latex','fontsize',15);

exportgraphics(gcf,'./gas_data/disp_field_step40.pdf','ContentType','vector');

% Plot stress with a colormap
% MPa
t1 = cellstress{plotstep}(1,:);
t2 = cellstress{plotstep}(2,:);
t3 = cellstress{plotstep}(4,:);
t4 = cellstress{plotstep}(3,:);

FontSize=16;
hh=figure('Name','Stress');
set(hh,'color','w');
% set(hh,'Position',[400 100 1000 800]);
subplot(2,2,1);
snscontour(node, element,fac, u_xp,u_yp,t1','{\it\sigma''}_{xx}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on


subplot(2,2,2);
snscontour(node, element,fac, u_xp,u_yp,t2','{\it\sigma''}_{yy}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on

subplot(2,2,3);
snscontour(node, element,fac, u_xp,u_yp,t3','{\it\sigma''}_{xy}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on

subplot(2,2,4);
snscontour(node, element,fac, u_xp,u_yp,t4','{\it\sigma''}_{zz}');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
hold on
sgtitle('Terzaghi effective stress {\boldmath$\sigma$''} components, MPa','Interpreter','latex','fontsize',15);

exportgraphics(gcf,'./gas_data/stress_field_step40.pdf','ContentType','vector');


% Plot equivalent plastic strain with a colormap
fac = 1;
figure
set(gcf,'Color','w');
t5=cellSDV{plotstep}(2,:);
snscontour(node,element,fac, u_xp,u_yp,t5','PEEQ');
c = colorbar; c.FontSize = 14; set(c,'TickLabelInterpreter','latex');
title('Equivalent plastic strain','Interpreter','latex','fontsize',15);
hold on

exportgraphics(gcf,'./gas_data/peeq_field_step40.pdf','ContentType','vector');


% Pore pressure
figure
set(gcf,'Color','w');
snscontour(node,element,fac,u_xp,u_yp,p_nodal,'PP');
c = colorbar; c.FontSize = 14; set(c,'TickLabelInterpreter','latex');
title('Gas pressure $p$, MPa','Interpreter','latex','fontsize',15);
hold on

exportgraphics(gcf,'./gas_data/gasp_field_step40.pdf','ContentType','vector');


p_at_corner = zeros(length(watch.dt),2);
time = 0;
% Gas pressure at corner
for step = 1:length(watch.dt)
    time = time + watch.dt(step);
    p_at_corner(step,:) = [time/24/3600, cellUP{step}(nnode*ndof+1)]; % (-1, -1) POINT; (1, 0) POINT: 135
end

figure;
set(gcf,'Color','w');
semilogx(p_at_corner(:,1),p_at_corner(:,2), 'b', 'Linewidth', 1); grid on; hold on;
% ref = load('./gas_data/reference_result.mat');
% semilogx(ref.COMSOL_ELASTIC(:,1),ref.COMSOL_ELASTIC(:,2), 'r-.', 'Linewidth', 1);
xlim([0.1, 500]);
xlabel('Time (d)', 'Interpreter','latex');
ylabel('Gas $p$ (MPa)', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca, 'FontSize', 13);
set(get(gca, 'XLabel'),'FontSize', 15);
set(get(gca, 'YLabel'),'FontSize', 15);

toc;

% save ./gas_data/reference_result.mat COMSOL_ELASTIC p_at_corner;

%% Compare with the case of incompressible fluid, pore pressure steady-state distribution (log)
elemType = 'T3';
poreliquid_ss = zeros(nnode,1);
pe = cellUP{plotstep}(nnode*ndof+1); % MPa, same as the corner gas pressure at Step 40 (66.6 days)
pw = 2; % MPa, BHP
Re = sqrt(2); % meter
Rw = 0.1; % meter, well radius
for i = 1:nnode
    r = sqrt(gcoord(i,1)^2 + gcoord(i,2)^2);
    poreliquid_ss(i) = pw + (pe - pw)/log(Re/Rw)*log(r/Rw);
end

figure;
plot_field(gcoord,ele_nods,elemType,poreliquid_ss);
colormap(jet);
c = colorbar; c.FontSize = 14; set(c,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',15);
title('Incompressible fluid steady-state $p$, MPa','Interpreter','latex','fontsize',15);

%% user-defined mesh function
function [hfun] = hfun8(test)
%HFUN8 user-defined mesh-size function for DEMO-8.

hmax = 0.1;
hmin = 0.02;

xmid = 0.0;
ymid = 0.0;

hcir = exp( -1*(test(:,1)-xmid).^2 ...
    -1*(test(:,2)-ymid).^2 );

hfun = hmax - (hmax-hmin) * hcir;

end
