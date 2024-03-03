% A mixed formulation of compressible elasticity capable of representing
% the incompressible limit

% Reference: The Finite Element Method - Thomas J. R. Hughes

% 单元体测试结果表明对于泊松比显著小于0.5的情况 精度不行 可能算法还是有缺陷
% 尝试sigma_ij = K*theta*delta_ij + 2*G*e_(i,j)
% 原先sigma_ij = lambda*theta*delta_ij + 2*G*e_ij
% 改Ce = 2*G*(diag([1,1,1,0.5,0.5,0.5])-1/3*(one_voigt*one_voigt')) 共两处 当前file和UMAT文件
% 改P_MEAN的计算方法
% 改data_const.mobility

restoredefaultpath;
clear; clearvars -global; clc; close all;
global ele_nods gcoord nnode nnel nel ndof ndofp bcdof bcval_incr g_const
format short;
addpath('./SFEM_basic/');
addpath(genpath('./dengwirda-mesh2d-ceb68eb'));
addpath('./plottingN');
addpath('./meshingN/');

%% Geometry and material of Cook's membrane
geoheight=10;  % meter
geowidth=10;
foot=1;

keypoint = [0, 0;
    geowidth, 0;
    geowidth, geoheight;
    0, geoheight;
    foot, geoheight;
    foot, 0;
    0, 0];

meshsize = [2, 0.04, 1, 9];

[gcoord, ele_nods] = meshfooting('Non-uniform mesh', keypoint, meshsize, true);

TRR = ST3Element(ele_nods, gcoord);
TRR.PlotElements;
hold on;

nnode = size(gcoord, 1);
nnel = 3;
ndof = 2;
ndofp = 1;
nel = size(ele_nods, 1);
g_const = 0;       % set it to zero could eliminate the effect of gravity

% kPa
E = 1000; nu = 0.499; tau = 2; rho_buo = 0; % buoyant density of the mixture
data_const.tauG = tau/(E/(1+nu));   % Stabilization multiplier, already divided by 2G
%data_const.mobility = (1+nu)*(1-2*nu)/E/nu; % 1/Lamé, Set delta_t = 1
data_const.mobility = 3*(1-2*nu)/E;

phi = 0; psi = 0; coh = 1; % cohesion in kPa
Props=ones(nnode, 1)*[E, nu, phi, psi, coh];

in_situ_stress = [0; 0; 0; 0; 0; 0]; % initial stress field, in general, the first three components could be negative
data_const.Ce_t0 = E/2/(1+nu)*diag([2,2,1]); % only contains the "mu" part

old_solution = sparse((ndof + ndofp)*nnode, 1);
new_solution = sparse((ndof + ndofp)*nnode, 1);
residual_traction = sparse(ndof*nnode, 1);
r_trac_old = residual_traction;

b_left=find(gcoord(:,1)==0);
b_bottom=find(gcoord(:,2)==0);
b_right=find(gcoord(:,1)==geowidth);
b_top=find(gcoord(:,2)==geoheight);
b_foot =intersect(find(gcoord(:,1)<=foot), b_top);

% Dirichlet BC
bcdof = [ndof*b_foot', ndof*b_foot'-1, ndof*b_left'-1, ndof*b_right'-1, ndof*b_bottom'-1, ndof*b_bottom']; % 位移加载
%bcdof = [ndof*b_left'-1, ndof*b_right'-1, ndof*b_bottom'-1, ndof*b_bottom'];
bcdof = unique(bcdof, 'stable');


%% 位移加载
Total_U = -0.05;
num_steps = 100;
watch.dt = 10*ones(1,num_steps);   % time STEP INTERVAL (pseudo time here, 10 is arbitrary!)
watch.now = 0; % start
bcval_incr = sparse(length(watch.dt), length(bcdof));
bcval_incr(1:num_steps, 1:length(b_foot)) = Total_U/num_steps;
epsp = 0.001;

traction_f = @(x,t)([0;0]); % no prescribed traction

%% 应力加载
% num_steps = 30;
% watch.dt = 1*ones(1,num_steps);   % time STEP INTERVAL (pseudo time here, 1 is arbitrary!)
% watch.now = 0; % start
% bcval_incr = sparse(length(watch.dt), length(bcdof));
% epsp = 0.01;
% max_t = sum(watch.dt);
% traction_f = @(x,t)([0;-6./max_t.*t.*(x<=foot)]); % kPa

%% FEM Solving (changing part)
% Find adjacent elements of each node
[nod_adjele] = get_nod_adjele;
% Compute the area of SD associated with node and element areas
[area_nod, area_T3] = cal_area_nod_T3(nod_adjele);

[K_example, all_sd_set_node, all_sd_B, all_sd_E]...
    = pre_assemble_BigK(nod_adjele, area_nod, area_T3);

%% all_sd_N: Integration of shape function on SD divided by SD area
% First assume psi = 0, thus all_sd_N would not be used
all_sd_N = cal_smooth_N(gcoord, ele_nods, all_sd_set_node, nod_adjele, area_nod, area_T3);

%% Continue
[Mass_mat, S_PPP, int_Nu] = pre_assemble_MassN(area_T3);
[K11_correction, ~] = assemble_stab(nod_adjele, area_nod, area_T3,...
    all_sd_set_node, all_sd_B, all_sd_E,...
    data_const, epsp);

sdof = (ndof + ndofp)*nnode;
index_not_constrained = setdiff(1:sdof, bcdof);
max_newton_iter = 300; tol = 5e-4; tol2 = 1e-7; % tol2 not used

cellstress = cell(length(watch.dt), 1); % Store stress tensor on every node for all time steps
cellUP = cell(length(watch.dt), 1); % Store the "new_solution" for all time steps
cellSDV = cell(length(watch.dt), 1);
stress = in_situ_stress*ones(1,nnode);  % Old
stress_new = stress;
stress_elastic = stress;
SDV = zeros(1,nnode); % Solution-Dependent State Variables
SDV_new = SDV;

%% For this problem with small strain assumption, K matrix would not change
K = assemble_system_K(area_nod, all_sd_set_node, all_sd_B, ...
    data_const, K_example, S_PPP, Mass_mat, Props);
% NS-FEM to SNS-FEM
K(1:nnode*ndof, 1:nnode*ndof) = K(1:nnode*ndof, 1:nnode*ndof) + K11_correction;

%RHS_correction = sparse(ndof*nnode, 1); % stablizied nodal integration over SD?

%% Time step loop
for step = 1:length(watch.dt)
    watch.now = watch.now + watch.dt(step);
    delta_t  = watch.dt(step);
    print_info = 'LOAD STEP = %d; TIME = %.2f:\n';
    fprintf(print_info, step, watch.now);
    r_trac_new = assign_tractionBC2(residual_traction, b_top, traction_f, watch.now);
    
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
        [RHS, stress_new, stress_elastic, SDV_new] = update_stress_ep(area_nod, ...
            all_sd_set_node, all_sd_B, all_sd_N, new_sol_incr, stress_new, SDV_new, Props, new_solution); % SDV and stress will change during iteration!

        %RHS_Urow = test_residual_Urow(area_nod, all_sd_set_node, all_sd_B, K_example, stress_elastic, new_solution);
        %RHS_Urow = RHS_Urow - K11_correction*new_solution(1:ndof*nnode);
        %ind_U_not_fix = setdiff(1:nnode*ndof, bcdof);
        %disp(norm(RHS_Urow(ind_U_not_fix))); % true residual of the U-row (balance of linear momentum) in mixed FEM

        print_info = '\t NEWTON = %d; RES = %.4E; \n';
        fprintf(print_info, iter, norm(RHS));
        
        if norm(RHS) < tol
            stress_new = stress_elastic; 
            % After Step 4, old "stress_new" is updated to stress_elastic, and pull back to yield surface
            % That is the new "stress_new" at line 163 LHS
            % This stress_elastic satisfies the momentum balance equation accurately!
            break;
        end

        % Solve and update
        new_solution_temp = new_solution; % Store the temporary value BEFORE Newton update
        new_solution(index_not_constrained) = new_solution(index_not_constrained)...
        + K(index_not_constrained, index_not_constrained)...
        \RHS(index_not_constrained);
        new_sol_incr = new_solution - new_solution_temp;
        error_temp = norm(new_sol_incr(1:nnode*ndof))/(norm(new_sol_incr0(1:nnode*ndof)) + eps);

        disp(['               Solve, U-err = ', num2str(error_temp)]);
        
        if iter == max_newton_iter
            error('Reach the largest NR iteration numbers!'); % Only first (step-1) results should be used!
        end
    end
    
    % Save state and store (converged) result
    %RHS_correction = K11_correction*(new_solution(1:ndof*nnode) - old_solution(1:ndof*nnode)); % for next time step

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

% save('../../footing_phi0psi0_vlocking_100steps.mat',"cellstress","cellSDV","cellUP","resi_force","Total_U");

%% Plot sigma or mean stress of Cook's membrane
u_xp = new_solution(1:2:nnode*ndof-1);
u_yp = new_solution(2:2:nnode*ndof);
p_nodal = new_solution(nnode*ndof+1:end); % [Pa]
t1 = stress(1,:);
t2 = stress(2,:); % deviatoric stress (Pa) yy component
t3 = stress(3,:);
t4 = stress(4,:); % tau_xy

node = gcoord;
element = ele_nods;

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, t1'-p_nodal, 'name');
h = colorbar;
t = get(h,'Limits');
h.Label.String = 'Total stress xx (kPa)';
axis on; box on;
clim([-4, 0.2]);
%clim([-8, 0.4]);
set(gca,'Fontsize',14,'FontName','Times new Roman');
colormap('parula');

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, t2'-p_nodal, 'name');
h = colorbar;
t = get(h,'Limits');
h.Label.String = 'Total stress yy (kPa)';
axis on; box on;
clim([-6, 0.5]);
%clim([-12, 0.5]);
set(gca,'Fontsize',14,'FontName','Times new Roman');

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, t3'-p_nodal, 'name');
h = colorbar;
%t = get(h,'Limits');
h.Label.String = 'Total stress zz (kPa)';
axis on; box on;
clim([-4, 0.5]);
%clim([-10, 0.5]);
set(gca,'Fontsize',14,'FontName','Times new Roman');

figure;
set(gcf, 'Color', 'w');
snscontour(node, element, 0, u_xp, u_yp, t4', 'name');
h = colorbar;
%t = get(h,'Limits');
h.Label.String = 'Shear stress xy (kPa)';
axis on; box on;
clim([-0.8, 1]);
set(gca,'Fontsize',14,'FontName','Times new Roman');


figure;
set(gcf, 'Color', 'w');
t5=cellSDV{end};
snscontour(node, element, 1, u_xp, u_yp, t5','PEEQ');
axis on; box on;
clim([0, 0.2]);
set(gca,'Fontsize',14,'FontName','Times new Roman'); hold on;

% plot_mesh(node, element, 'T3', 'k-');

%% Results not good for stress: stress singularity at the footing edge!


%% Plot normalized vertical reaction force versus penetration depth
resi_force = zeros(1,num_steps);
temp = [gcoord(b_foot, 1), b_foot];
temp = sortrows(temp);
b_foot = temp(:, 2); % order nodes in b_foot by x-coord
clearvars temp
for step = 1:num_steps
    resi_force(step) = trapz(gcoord(b_foot, 1),cellstress{step}(2, b_foot')'-cellUP{step}(nnode*ndof+b_foot));
end
figure;
set(gcf, 'Color', 'w');
plot(-Total_U/num_steps*(0:1:num_steps), [0,-resi_force]/coh/foot/1, 'b-d','Linewidth', 1, 'MarkerFaceColor','y'); grid on; % thickness in z direction: 1 unit
xlabel('Penetration depth (m)','Interpreter', 'LaTeX');
ylabel('Normalized Reaction force $q/c_{\rm u}$', 'Interpreter', 'LaTeX');
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


%% User-defined mesh function for plate with a circular hole
function [hfun] = hfun8(test)

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

function [RHS, stress_new, stress_elastic, SDV_new] = update_stress_ep(area_nod, ...
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
    [stress_new(:, ino), stress_elastic(:,ino), SDV_new(:, ino), epvol_incr, stress_diff] = ...
        DP_modified(Props(ino,:), stress(:, ino), strain_ino_incr, SDV(:, ino), new_sol(nnode*ndof+ino)); % last input for mixed pressure
    
    % stress_new is inside or on the yield surface
    residual_1(index_u) = residual_1(index_u) +...
        transpose(all_sd_B{ino})*stress_diff([1,2,4])*area_nod(ino);
    residual_2(index_p) = residual_2(index_p) + transpose(all_sd_N{ino})*epvol_incr*area_nod(ino);

    RHS = [residual_1; -residual_2];
end

end

function [STRESS, STRESS_INTERM, hsv, epvol_incr, STRESS_DIFF] = DP_modified(PROPS, STRESS0, DSTRAIN0, hsv0, P_MEAN)
% Input P_MEAN is + for compression

DSTRAIN0 = [DSTRAIN0(1); DSTRAIN0(2); 0; DSTRAIN0(3); 0; 0]; % 6*1 VECTOR

E = PROPS(1); %  Young's modulus
nu = PROPS(2); %  Poisson's ratio
G = E/2/(1+nu);
phi = PROPS(3);  %  phi  -  angle of friction (rad)
psi = PROPS(4);  %  psi  -  angle of dilation (rad)
coh = PROPS(5);  %  c  -  cohesion

one_voigt = [1,1,1,0,0,0]';
Ce = 2*G*(diag([1,1,1,0.5,0.5,0.5])-1/3*(one_voigt*one_voigt')); % change between 0 and 1

A = 3*sqrt(2)*coh*cos(phi)/sqrt(9+3*(sin(phi))^2);
B = 3*sqrt(2)*sin(phi)/sqrt(9+3*(sin(phi))^2);
b = 3*sqrt(2)*sin(psi)/sqrt(9+3*(sin(psi))^2);

STRESS_INTERM = STRESS0 + Ce*DSTRAIN0; % Intermediate deviatoric stress quantity
stress_trial_tensor = [STRESS_INTERM(1), STRESS_INTERM(4), STRESS_INTERM(5);
    STRESS_INTERM(4), STRESS_INTERM(2), STRESS_INTERM(6);
    STRESS_INTERM(5), STRESS_INTERM(6), STRESS_INTERM(3)];

stress_dev_trial = stress_trial_tensor - 1/3*trace(stress_trial_tensor)*eye(3);
q_trial = sqrt(1.5)*norm(stress_dev_trial, 'fro');

%P_MEAN = -P_MEAN*(3*nu/(1+nu)); % change sign convention; convert to TRUE mean stress
P_MEAN = -P_MEAN;
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

function mat_N_SD = cal_smooth_N(nodeT3, elementT3, set_nod_SD, adjele_nod, area_SD, area_T3)

nnodeT3 = height(nodeT3);
mat_N_SD = cell(nnodeT3,1);
for ino = 1:nnodeT3
	adjele_ino = adjele_nod{ino};
	nodset_ino = set_nod_SD{ino};
	N_ino = zeros(1,length(nodset_ino));
	for iel = 1:length(adjele_ino)
		nodeT3_iel = elementT3(adjele_ino(iel),:);
		n1 = nodeT3_iel(1);
		n2 = nodeT3_iel(2);
		n3 = nodeT3_iel(3);
		
		x1 = nodeT3(n1,1);   y1 = nodeT3(n1,2);
		x2 = nodeT3(n2,1);   y2 = nodeT3(n2,2);
		x3 = nodeT3(n3,1);   y3 = nodeT3(n3,2);

        xm1 = (x1+x2)/2;      ym1 = (y1+y2)/2;
        xm2 = (x2+x3)/2;      ym2 = (y2+y3)/2;
        xm3 = (x3+x1)/2;      ym3 = (y3+y1)/2;
        xc = (x1+x2+x3)/3;    yc = (y1+y2+y3)/3;
		
        if ino == n1
            xc_sub = (x1+xm1+xc+xm3)/4;
            yc_sub = (y1+ym1+yc+ym3)/4;
        elseif ino == n2
            xc_sub = (x2+xm2+xc+xm1)/4;
            yc_sub = (y2+ym2+yc+ym1)/4;            
        elseif ino == n3
            xc_sub = (x3+xm3+xc+xm2)/4;
            yc_sub = (y3+ym3+yc+ym2)/4;            
        end
        
		% calculate factor of shape function
        a = 1/(2*area_T3(adjele_ino(iel)))*[x2*y3-x3*y2, x3*y1-x1*y3, x1*y2-x2*y1];   % a = [a1,a2,a3]
        b = 1/(2*area_T3(adjele_ino(iel)))*[y2-y3, y3-y1, y1-y2];                     % b = [b1,b2,b3] d/dx
        c = 1/(2*area_T3(adjele_ino(iel)))*[x3-x2, x1-x3, x2-x1];                     % c = [c1,c2,c3] d/dy
		
		NT3_subc = a + b*xc_sub + c*yc_sub;
		N_gp = (NT3_subc)*area_T3(adjele_ino(iel))/3/area_SD(ino);
		
        % assemble N_ino
        for iino = 1:length(nodeT3_iel) % iino = 1:3
            iinode = nodeT3_iel(iino);
            seq = find(nodset_ino==iinode);
            seq = seq(1);
            N_ino(1,seq) = N_ino(1,seq) + N_gp(iino);
        end
	end % end of loop over adjacent elements
	mat_N_SD{ino} = N_ino;
end

end

%% For testing some convergence behaviors
function RHS_Urow = test_residual_Urow(area_nod, all_sd_set_node, all_sd_B, K_example, stress_new, new_solution)

global ndof nnode
residual_1 = sparse(ndof*nnode, 1);

for ino = 1:nnode
    nodB = all_sd_set_node{ino};
    index_u = get_eledof(nodB, length(nodB), ndof);
    residual_1(index_u) = residual_1(index_u) - transpose(all_sd_B{ino})*stress_new([1,2,4], ino)*area_nod(ino);
end
K12 = K_example(1:ndof*nnode, ndof*nnode+1:end);
RHS_Urow = residual_1 + K12*new_solution(ndof*nnode+1:end);

end
