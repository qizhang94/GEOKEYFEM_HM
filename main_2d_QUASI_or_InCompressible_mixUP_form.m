% A mixed formulation of compressible elasticity capable of representing
% the incompressible limit

% Reference: The Finite Element Method - Thomas J. R. Hughes

restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel
global ndof ndofp bcdof bcval_incr g_const
format short;
addpath('./SFEM_basic/');
addpath(genpath('./dengwirda-mesh2d-ceb68eb'));
addpath('./plottingN');
addpath('./meshingN/');

%% Geometry of Infinite plate with a circular hole
% node = [0, 5;
%     5, 5;
%     5, 0];
% edge = [0, 1;
%     1, 2;
%     2, 3;
%     3, 1];
% 
% adel = 2.*pi / +32;
% amin = 0;
% amax = pi/2.;
% 
% xcir = +1. * ...
%     cos(amin:adel:amax)';
% ycir = +1. * ...
%     sin(amin:adel:amax)';
% ncir = [xcir,ycir];
% ncir(end,1) = 0;
% numc = size(ncir,1); % number of keypoints in the inner circle
% 
% ecir(:,1) = ...
%     (1:numc-1)';
% ecir(:,2) = ...
%     (2:numc-0)';
% 
% node = [ncir; node];
% 
% edge = [ecir; edge+numc];
% edge(end,2) = 1;
% 
% %---------------------------------------------- do mesh-gen.
% hfun = @hfun8;
% 
% [vert,etri, ...
%     tria,tnum] = refine2(node,edge,[],[],hfun);
% 
% %---------------------------------------------- do mesh-opt.
% [vert,~, ...
%     tria,~] = smooth2(vert,etri,tria,tnum);
% 
% gcoord = vert; ele_nods = tria;
% 
% TRR = ST3Element(ele_nods, gcoord);
% TRR.PlotElements;
% hold on;

%% Material properties, boundary condition of Infinite plate with a circular hole 
% nnode = size(gcoord, 1);
% nnel = 3;
% ndof = 2;
% ndofp = 1;
% nel = size(ele_nods, 1);
% g_const = 0;       % set it to zero could eliminate the effect of gravity
% 
% % Pa
% E = 1000; nu = 0.4999999; tau = 2; rho_buo = 0; % buoyant density of the mixture
% data_const.tauG = tau/(E/(1+nu));   % Stabilization multiplier, already divided by 2G
% data_const.mobility = (1+nu)*(1-2*nu)/E/nu; % Set delta_t = 1
% 
% Props=ones(nnode, 1)*[E, nu];
% 
% in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
% [~, ~, cto_ela_t0] =...
%         LinEla_mupart_UMAT(Props(1,:), in_situ_stress, zeros(3,1), 0);
% data_const.Ce_t0 = cto_ela_t0; % only contains the "mu" part
% 
% 
% old_solution = sparse((ndof + ndofp)*nnode, 1);
% new_solution = sparse((ndof + ndofp)*nnode, 1);
% residual_traction = sparse(ndof*nnode, 1);
% 
% b_left=find(gcoord(:,1)==0);
% b_bottom=find(gcoord(:,2)==0);
% b_right=find(gcoord(:,1)==5);
% 
% % Dirichlet BC
% bcdof = [ndof*b_left'-1, ndof*b_bottom'];
% bcdof = unique(bcdof, 'stable');
% 
% % Neumann BC
% 
% watch.dt = ones(1,1);   % Only one step
% watch.now = 0; % start
% bcval_incr = sparse(length(watch.dt), length(bcdof));
% epsp = 1;
% 
% traction_f = @(x,t)([1*t/sum(watch.dt);0]); % Tension on the right 1[Pa]

%% Geometry and material of Cook's membrane
numx = 16;
numy = 16;
[gcoord,ele_nods] = mesh_region([0,0], [48,44], [48,60], [0,44], numx,numy,'T3');

% TRR = ST3Element(ele_nods, gcoord);
% TRR.PlotElements;
% hold on;

nnode = size(gcoord, 1);
nnel = 3;
ndof = 2;
ndofp = 1;
nel = size(ele_nods, 1);
g_const = 0;       % set it to zero could eliminate the effect of gravity

% Pa
E = 1; nu = 0.4999999; tau = 2; rho_buo = 0; % buoyant density of the mixture
data_const.tauG = tau/(E/(1+nu));   % Stabilization multiplier, already divided by 2G
data_const.mobility = (1+nu)*(1-2*nu)/E/nu; % 1/Lamé, Set delta_t = 1

Props=ones(nnode, 1)*[E, nu];

in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
[~, ~, cto_ela_t0] =...
        LinEla_mupart_UMAT(Props(1,:), in_situ_stress, zeros(3,1), 0);
data_const.Ce_t0 = cto_ela_t0; % only contains the "mu" part


old_solution = sparse((ndof + ndofp)*nnode, 1);
new_solution = sparse((ndof + ndofp)*nnode, 1);
residual_traction = sparse(ndof*nnode, 1);

b_left=find(gcoord(:,1)==0);
b_right=find(gcoord(:,1)==48);

% Dirichlet BC
bcdof = [ndof*b_left'-1, ndof*b_left']; % clampled left edge
bcdof = unique(bcdof, 'stable');

% Neumann BC

watch.dt = ones(1,1);   % Only one step
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
epsp = 1;

traction_f = @(x,t)([0;1/16]); % Shear on the right 1[Pa]

%% FEM Solving
% Find adjacent elements of each node
[nod_adjele] = get_nod_adjele;
% Compute the area of SD associated with node and element areas
[area_nod, area_T3] = cal_area_nod_T3(nod_adjele);

[K_example, all_sd_set_node, all_sd_B, all_sd_E] = pre_assemble_BigK(nod_adjele, area_nod, area_T3);
[Mass_mat, S_PPP, int_Nu] = pre_assemble_MassN(area_T3);
[K11_correction, ~] = assemble_stab(nod_adjele, area_nod, area_T3,...
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
    residual_traction = assign_tractionBC2(residual_traction, b_right, traction_f, watch.now, 'vertical');

    for iter = 1:max_newton_iter

        
        [K, RHS, stress_new, SDV_new] = assemble_system(nod_adjele, area_nod, area_T3,...
        all_sd_set_node, all_sd_B, all_sd_E, ...
        new_solution, old_solution,...
        data_const, stress, SDV,...
        K_example, S_PPP, Mass_mat, int_Nu, residual_traction, delta_t, ...
        Props, true); % NS-FEM

        % NS-FEM to SNS-FEM
        K(1:nnode*ndof, 1:nnode*ndof) = K(1:nnode*ndof, 1:nnode*ndof) + K11_correction;
        RHS(1:nnode*ndof) = RHS(1:nnode*ndof) - K11_correction*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode));

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
    
end

clearvars K_example Mass_mat S_PPP int_Nu
clearvars ino step print_info delta_t error_temp sdof
clearvars iter
clearvars new_solution_temp

%% Plot sigma_yy of Infinite plate with a circular hole 
% u_xp = new_solution(1:2:nnode*ndof-1);
% u_yp = new_solution(2:2:nnode*ndof);
% p_nodal = new_solution(nnode*ndof+1:end);
% t2 = stress(2,:); % deviatoric stress
% 
% node = gcoord;
% element = ele_nods;
% 
% figure;
% snscontour(node, element, 0, u_xp,u_yp, t2'-p_nodal,'{\it\sigma}_{yy}');
% clim([-1,0.5]);
% set(gca,'Fontsize',14,'FontName','Times new Roman');

%% Plot sigma_xx of Cook's membrane
u_xp = new_solution(1:2:nnode*ndof-1);
u_yp = new_solution(2:2:nnode*ndof);
p_nodal = new_solution(nnode*ndof+1:end); % [Pa]
t1 = stress(1,:); % deviatoric stress (Pa)

node = gcoord;
element = ele_nods;

figure;
snscontour(node, element, 0, u_xp,u_yp, -p_nodal, 'Mean stress');
h = colorbar;
t = get(h,'Limits');
clim([-0.28,0.11]);
set(gca,'Fontsize',14,'FontName','Times new Roman');

% The deflections of center of right side of the Cook's membrane problem
findnod = intersect(b_right, find(gcoord(:,2) == 52));
disp(new_solution(ndof*findnod));


%% user-defined mesh function
function [hfun] = hfun8(test)
%HFUN8 user-defined mesh-size function for DEMO-8.

hmax = 0.2;
hmin = 0.05;

xmid = 0.0;
ymid = 0.0;

hcir = exp( -0.1*(test(:,1)-xmid).^2 ...
    -0.1*(test(:,2)-ymid).^2 );

hfun = hmax - (hmax-hmin) * hcir;

end

