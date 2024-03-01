% A mixed formulation of compressible elasticity capable of representing
% the incompressible limit

% Reference: The Finite Element Method - Thomas J. R. Hughes

restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel ndof ndofp bcdof bcval_incr g_const
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
numx = 32;
numy = 32;
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

phi = pi/12; psi = 0; coh = 0.1; 
Props=ones(nnode, 1)*[E, nu, phi, psi, coh];

in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
data_const.Ce_t0 = E/2/(1+nu)*diag([2,2,1]); % only contains the "mu" part


old_solution = sparse((ndof + ndofp)*nnode, 1);
new_solution = sparse((ndof + ndofp)*nnode, 1);
residual_traction = sparse(ndof*nnode, 1);
r_trac_old = residual_traction;

b_left=find(gcoord(:,1)==0);
b_right=find(gcoord(:,1)==48);

% Dirichlet BC
bcdof = [ndof*b_left'-1, ndof*b_left']; % clampled left edge
bcdof = unique(bcdof, 'stable');

% Neumann BC

watch.dt = [0.5, 0.5];   % Two time steps
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
epsp = 0.1;

traction_f = @(x,t)([0;1/16*t]); % Shear on the right +1[N]

%% FEM Solving (changing part)
% Find adjacent elements of each node
[nod_adjele] = get_nod_adjele;
% Compute the area of SD associated with node and element areas
[area_nod, area_T3] = cal_area_nod_T3(nod_adjele);

[K_example, all_sd_set_node, all_sd_B, all_sd_E]...
    = pre_assemble_BigK(nod_adjele, area_nod, area_T3); % all_sd_N: Integration of shape function on SD
%% all_sd_N
% First assume psi = 0, thus all_sd_N would not be used
all_sd_N = cell(nnode, 1);
for ino = 1:nnode
    all_sd_N{ino} = zeros(1,length(all_sd_set_node{ino}));
end

%% Continue
[Mass_mat, S_PPP, int_Nu] = pre_assemble_MassN(area_T3);
[K11_correction, ~] = assemble_stab(nod_adjele, area_nod, area_T3,...
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

%% For this problem with small strain assumption, K matrix would not change
K = assemble_system_K(area_nod, all_sd_set_node, all_sd_B, ...
    data_const, K_example, S_PPP, Mass_mat, Props);
% NS-FEM to SNS-FEM
K(1:nnode*ndof, 1:nnode*ndof) = K(1:nnode*ndof, 1:nnode*ndof) + K11_correction;

%% Time step loop
for step = 1:length(watch.dt)
    watch.now = watch.now + watch.dt(step);
    delta_t  = watch.dt(step);
    print_info = 'LOAD STEP = %d; TIME = %.2f:\n';
    fprintf(print_info, step, watch.now);
    r_trac_new = assign_tractionBC2(residual_traction, b_right, traction_f, watch.now, 'vertical');
    
    % Suppose we get the residual_traction_incr (We need the traction increment now!)
    %% Elastic trial FEM calculation (not on the constitutive level)
    RHS = [r_trac_new-r_trac_old; sparse(ndofp*nnode, 1)];
    % Adjust RHS according to bcval_incr, ONLY has an impact for non-zero Dirichlet BC
    RHS = RHS - K(:, bcdof)*bcval_incr(step,:)';
    new_solution(bcdof) = old_solution(bcdof) + bcval_incr(step,:)';
    new_solution(index_not_constrained) = new_solution(index_not_constrained)...
        + K(index_not_constrained, index_not_constrained)...
        \RHS(index_not_constrained);
    new_sol_incr0 = new_solution - old_solution;
    new_sol_incr = new_sol_incr0;
    %% Plastic correction FEM calculation
    for iter = 1:max_newton_iter
        % 1. Calculate current p and deviatoric stress
        % 2. Return mapping to the yield surface, and assemble the
        % deviatoric stress differnce between trial quantity and corrected quantity, check the
        % convergence criterion
        % 3. Assemble the volumetric plastic strain for the second equation
        % 4. Solve and update "new_solution", and back to 1.
        [RHS, stress_new, SDV_new] = update_stress_ep(area_nod, ...
            all_sd_set_node, all_sd_B, all_sd_N, new_sol_incr, stress_new, SDV_new, Props, new_solution); % SDV and stress will change during iteration!
        
        print_info = '\t NEWTON = %d; RES = %.4E; \n';
        fprintf(print_info, iter, norm(RHS));
        
        if norm(RHS) < tol
            break;
        end

        % Solve and update
        new_solution_temp = new_solution; % Store the temporary value BEFORE Newton update
        new_solution(index_not_constrained) = new_solution(index_not_constrained)...
        + K(index_not_constrained, index_not_constrained)...
        \RHS(index_not_constrained);
        new_sol_incr = new_solution - new_solution_temp;
        error_temp = norm(new_sol_incr(1:nnode*ndof))/(norm(new_sol_incr0(1:nnode*ndof)) + eps);

        disp(['                  U-err = ', num2str(error_temp)]);
        
    end
    
    % Save state and store (converged) result
    stress = stress_new;
    old_solution = new_solution;
    SDV = SDV_new;
    cellstress{step} = stress;
    cellUP{step} = old_solution;
    cellSDV{step} = SDV;
    r_trac_old = r_trac_new;
    
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

%% Plot sigma or mean stress of Cook's membrane
u_xp = new_solution(1:2:nnode*ndof-1);
u_yp = new_solution(2:2:nnode*ndof);
p_nodal = new_solution(nnode*ndof+1:end); % [Pa]
t1 = stress(1,:); % deviatoric stress (Pa)

node = gcoord;
element = ele_nods;

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, t1'-p_nodal, 'name');
h = colorbar;
t = get(h,'Limits');
h.Label.String = 'Total stress xx';
axis on; box on;
%clim([-0.28,0.11]);
set(gca,'Fontsize',14,'FontName','Times new Roman');

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, -p_nodal, 'name');
h = colorbar;
t = get(h,'Limits');
h.Label.String = 'Mean stress';
axis on; box on;
%clim([-0.28,0.11]);
set(gca,'Fontsize',14,'FontName','Times new Roman');

figure;
set(gcf, 'Color', 'w');
t5=cellSDV{end};
snscontour(node, element, 1, u_xp, u_yp, t5','PEEQ');
axis on; box on;
set(gca,'Fontsize',14,'FontName','Times new Roman'); hold on;
plot_mesh(node, element, 'T3', 'k-');


% The deflections of center of right side of the Cook's membrane problem
findnod = intersect(b_right, find(gcoord(:,2) == 52));
disp(['Right side center deflection = ', num2str(new_solution(ndof*findnod)), '(m)']);


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

function K_UP = assemble_system_K(area_nod, all_sd_set_node, all_sd_B, ...
    material_data_const, K_example, S_PPP, Mass_mat, Props)

global ndof nnode

K11 = sparse(ndof*nnode, ndof*nnode);

for ino = 1:nnode
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    
    % Constitutive model **************************************************
    [~, ~, cto] = LinEla_mupart_UMAT(Props(ino,:), zeros(6,1), [0;0;0], 0);

    % Take care of the "sign" (+ or -)
    K11(index_u, index_u) = K11(index_u, index_u) + transpose(all_sd_B{ino})*cto*all_sd_B{ino}*area_nod(ino);
  
end
K12 = -K_example(1:ndof*nnode, ndof*nnode+1:end);   % K21 = transpose(K12) because of symmetry

factor = material_data_const.tauG;
K22_uns = -Mass_mat*(material_data_const.mobility); % For nearly incompressibility u-p form
K22 = K22_uns - factor*S_PPP;    % Pressure stabilization, factor = \tau/(2*G)

K_UP = [K11, K12; transpose(K12), K22]; % K22 semi-negative definite

end

function [RHS, stress_new, SDV_new] = update_stress_ep(area_nod, ...
    all_sd_set_node, all_sd_B, all_sd_N, new_sol_incr, stress, SDV, Props, new_sol) % _ep stands for elastoplasticity 

global ndof ndofp nnode

residual_1 = sparse(ndof*nnode, 1);
residual_2 = sparse(ndofp*nnode, 1);

for ino = 1:nnode
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    index_p = get_eledof(nodB, length(nodB), ndofp);
    strain_ino_incr = all_sd_B{ino}*new_sol_incr(index_u);

    % Modified Drucker-Prager model
    [stress_new(:, ino), SDV_new(:, ino), epvol_incr, stress_diff] = ...
        DP_modified(Props(ino,:), stress(:, ino), strain_ino_incr, SDV(:, ino), new_sol(nnode*ndof+ino)); % last input for mixed pressure
    
    % stress_new is inside or on the yield surface
    residual_1(index_u) = residual_1(index_u) +...
        transpose(all_sd_B{ino})*stress_diff([1,2,4])*area_nod(ino);
    residual_2(index_p) = residual_2(index_p) + transpose(all_sd_N{ino})*epvol_incr*area_nod(ino);

    RHS = [residual_1; -residual_2];
end

end

function [STRESS, hsv, epvol_incr, STRESS_DIFF] = DP_modified(PROPS, STRESS0, DSTRAIN0, hsv0, P_MEAN)
% Input P_MEAN is + for compression

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR

E = PROPS(1); %  Young's modulus
nu = PROPS(2); %  Poisson's ratio
G = E/2/(1+nu);
phi = PROPS(3);  %  phi  -  angle of friction (rad)
psi = PROPS(4);  %  psi  -  angle of dilation (rad)
coh = PROPS(5);  %  c  -  cohesion

Ce = G*diag([2,2,2,1,1,1]);

A = 3*sqrt(2)*coh*cos(phi)/sqrt(9+3*(sin(phi))^2);
B = 3*sqrt(2)*sin(phi)/sqrt(9+3*(sin(phi))^2);
b = 3*sqrt(2)*sin(psi)/sqrt(9+3*(sin(psi))^2);

STRESS_INTERM = STRESS0 + Ce*DSTRAIN0; % Intermediate deviatoric stress quantity
stress_trial_tensor = [STRESS_INTERM(1), STRESS_INTERM(4), STRESS_INTERM(5);
    STRESS_INTERM(4), STRESS_INTERM(2), STRESS_INTERM(6);
    STRESS_INTERM(5), STRESS_INTERM(6), STRESS_INTERM(3)];

stress_dev_trial = stress_trial_tensor - 1/3*trace(stress_trial_tensor)*eye(3);
q_trial = sqrt(1.5)*norm(stress_dev_trial, 'fro');

P_MEAN = -P_MEAN*(3*nu/(1+nu)); % change sign convention; convert to TRUE mean stress
if (sqrt(2/3)*q_trial + B*P_MEAN - A <= 0)
    hsv = hsv0;
    STRESS = STRESS_INTERM;
    epvol_incr = 0;
    STRESS_DIFF = zeros(6,1);
else
    nhat_tensor = stress_dev_trial/norm(stress_dev_trial, 'fro');
    nhat = [nhat_tensor(1,1); nhat_tensor(2,2); nhat_tensor(3,3);...
        nhat_tensor(1,2); nhat_tensor(1,3); nhat_tensor(2,3)];
    q_con = (A - B*P_MEAN)/sqrt(2/3);
    if q_con < 0
        q_con = 0;
    end
    STRESS_DIFF = sqrt(2/3)*(q_trial - q_con)*nhat;
    STRESS = STRESS_INTERM - STRESS_DIFF;
    dlambda = (q_trial - q_con)/sqrt(6)/G;
    hsv = hsv0 + sqrt(2/3)*dlambda;
    epvol_incr = dlambda*b;
end

end
