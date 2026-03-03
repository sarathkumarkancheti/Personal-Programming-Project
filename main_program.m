%% This MATLAB program solves the FEM formulation for the problems discussed
%% in Liu, Y., Zhang, H., Zhang, J. and Zheng, Y., 2016....
%% Transient swelling of polymeric hydrogels: A new finite element solution framework....
%% International Journal of Solids and Structures, 80, pp.246-260.[https://doi.org/10.1016/j.ijsolstr.2015.11.010]
%% __________________________________________________________________________
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% main program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nnode_ele_m              : No.of displacement nodes per element
% nnode_ele_c              : No.of chemical potential nodes per element
% nnodes_m                 : No.of displacement nodes in the mesh
% nnodes_c                 : No.of chemical potential nodes in the mesh
% nodeCoords               : Nodal coordinates (x,y) for 2D , (x,y,z) for 3D
% ncoord                   : No.of coordinates for each node 2 for 2D, 3 for 3D
% elementNodes             : Element nodal connectivity matrix
% flux_nodes               : Nodes on which chemical boundary conditions are applied
% numElements              : No.of elements in the mesh
% cp_solvent               : External solvent chemical potential
% j_flux                   : External fluid flux
% m_d                      : Nominal diffusion constant               
% dt                       : Time increment Δt
% tend                     : End of simulation time
% t0                       : Initial time
% nt                       : No.of time steps (t0-tend)/Δt
% cp0                      : Intial chemical potential [μ]_0
% max_iterations           : Maximum No.of iterations for Newton-Raphson loop
% stretch                  : Stretch of all elements for all time steps
% sigma_xx                 : Stress in xx direction at integration points
% sigma_yy                 : Stress in yy direction at integration points
% sigma_xy                 : Stress in xy direction at integration points
% U_global                 : DOF matrix [w;μ]
% U_nr                     : Newton iteration DOF matrix
% du                       : Incremental dof Δu
% error                    : Relative error for tolerance check |Δu_k|/|U_nr_k+1|
% element_displacement     : Nodal displacement of an element [w] = [ux;uy] for 2D , [ux;uy;uz] for 3D
% ele_incr                 : Incremental nodal displacement of an element [Δux;Δuy] for 2D, [Δux;Δuy;Δuz] for 3D
% element_cp               : Nodal chemical potential of an element [μ]
% boundary_local_nodes     : Boundary local node indices of an element on which external fluid flux or solvent is applied
% global_dofs              : displacement degree of freedom [ux;uy] for 2D , [ux;uy;uz] for 3D
% TestCase                 : To perform test cases make this 'Yes' and others as 'No' for any problem
% PatchTest                : To perform patch test make this 'Yes', problem = 1, TestCase = 'No', deformControl = 'No'
% deformControl            : To perform deformation control test make this 'Yes' , problem = 2 and others 'No'
%__________________________________________________________________________
clear;
problem = input('Enter the problem number: ');
%problem = 4;
TestCase = 'No';
patchTest = 'No';
deformControl = 'No';

%Initial normalized chemical potential calculated form F = I, stress = 0, X = 0.2
cp0 = (log(1-0.999)+1+0.2);

% Initialize model parameters
if problem == 1
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\convergence_1.txt', 'w');
ncoord = 2;
if strcmpi(patchTest, 'yes')
% Mesh data for patch test (problem = 6 in inp_coordinates function)
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                             right_edge_nodes,top_edge_nodes,nodeCoords,~] = inp_coordinates(2,6);
else
%Mesh data for 1D problem
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                             right_edge_nodes,top_edge_nodes,nodeCoords,~] = inp_coordinates(ncoord,problem);
end
nnode_ele_m = 4; 
nnode_ele_c = 4; 
nnodes_m = nnodes;
nnodes_c = nnodes;
numElements = length(elementIDs);
m_d = 0.5;
j_flux = 0; 
cp_solvent = 0;
dt = 0.005;
tend = 5;
nt = (tend-0)/dt+1;
stretch = ones(numElements,nt);
flux_nodes = nodeIDs(nodeCoords(:,2)==1);

% Initialize global degree of freedom
U_global = zeros(nnodes_m*ncoord+nnodes_c,nt);
% Initial chemical potential
U_global(nnodes_m*ncoord+1:end,1) = cp0; %[μ]_0

if strcmpi(TestCase, 'No') &&  strcmpi(patchTest, 'No')
% For TestCase = 'Yes' no need to prescribe chemical potential
U_global(nnodes_m*ncoord+flux_nodes,1) = cp_solvent; % cp_solvent on the boundary
stretch(20,1) = 1.498;
end

elseif problem == 2
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\convergence_2.txt', 'w');
ncoord = 2;
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                    right_edge_nodes,top_edge_nodes,nodeCoords,flux_nodes] = inp_coordinates(ncoord,problem);
nnode_ele_m = 8;
nnode_ele_c = 8;
nnodes_m = nnodes;
nnodes_c = nnodes;
numElements = length(elementIDs);
cp_solvent = 0;
j_flux = -4.3137e27;
m_d = 0.1;
dt = 0.005;
tend = 5;
nt = (tend-0)/dt+1;
stretch = ones(numElements,nt);
sigma_xx = zeros(numElements,4,nt);
sigma_yy = zeros(numElements,4,nt);
sigma_xy = zeros(numElements,4,nt);

% Initialize global degree of freedom
U_global = zeros(nnodes_m*ncoord+nnodes_c,nt);
% Initial chemical potential
U_global(nnodes_m*ncoord+1:end,1) = cp0; %[μ]_0

if strcmpi(deformControl,'Yes')
   % Perform simulation for one element
   nnode_ele_m = 8;
   nnode_ele_c = 8;
   nnodes_m = 8;
   nnodes_c = 8;
   numElements = 1;
   cp_solvent = 0;
   j_flux = 0; % No external flux for deformation control
   m_d = 0.1;
   dt = 0.005;
   elementNodes = [1,2,3,4,5,6,7,8];
   flux_nodes = [3,4,7];
   % Prescribe displacements
   Lambda = linspace(1,1.2,200);
   x_coords = [0;1;1;0;0.5;1;0.5;0];
   y_coords = [0;0;1;1;0;0.5;1;0.5];
   nodeCoords = [x_coords,y_coords];
   
   % Initialize global degree of freedom
   U_global = zeros(nnodes_m*ncoord+nnodes_c,nt);
   % Initial chemical potential
   U_global(nnodes_m*ncoord+1:end,1) = cp0; %[μ]_0
end

elseif problem == 3
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\convergence_3.txt', 'w');
ncoord = 3;
% Mesh data
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                          right_edge_nodes,top_edge_nodes,nodeCoords,~] = inp_coordinates(ncoord,problem);
nnode_ele_m = 8;
nnode_ele_c = 8;
nnodes_m = nnodes;
nnodes_c = nnodes;
numElements = length(elementIDs);
cp_solvent = 0;
j_flux = -1.1895e28;
m_d = 0.5;
dt = 0.005;
tend = 3;
nt = (tend-0)/dt+1;
flux_nodes = nodeIDs(nodeCoords(:,1)<=0.600000024 & nodeCoords(:,2)<=0.600000024 &nodeCoords(:,3)==1 ).';
xoz_fixed =  nodeIDs(nodeCoords(:,2)==0).';
yoz_fixed = nodeIDs(nodeCoords(:,1)==0).';
z_a = nodeIDs(nodeCoords(:,1)==0 & nodeCoords(:,2)==0 &nodeCoords(:,3)==1 );
z_b = nodeIDs(nodeCoords(:,1)==1 & nodeCoords(:,2)==1 &nodeCoords(:,3)==1 );
z_fixed_node = nodeIDs(nodeCoords(:,1)==1 & nodeCoords(:,2)==1 &nodeCoords(:,3)==0);

% Initialize global degree of freedom
U_global = zeros(nnodes_m*ncoord+nnodes_c,nt);
% Initial chemical potential
U_global(nnodes_m*ncoord+1:end,1) = cp0; %[μ]_0

elseif problem == 4
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\convergence_4.txt', 'w');
ncoord = 2;
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                            right_edge_nodes,top_edge_nodes,nodeCoords,~]= inp_coordinates(ncoord,problem);
nnode_ele_m = 8;
nnode_ele_c = 4;
nnodes_m = nnodes;
nnodes_c = 36;
flux_nodes = [top_edge_nodes(1:6),right_edge_nodes(1:6)];
numElements = length(elementIDs);
cp_solvent = 0;
j_flux = 0;
m_d = 0.01;
dt = 0.05;
tend = 200;
nt = (tend-0)/dt+1;
stretch = ones(numElements,nt);
sigma_xx = zeros(numElements,4,nt);
sigma_yy = zeros(numElements,4,nt);
sigma_xy = zeros(numElements,4,nt);
    
% Initialize global degree of freedom
U_global = zeros(nnodes_m*ncoord+nnodes_c,nt); 
% Initial chemical potential [μ]_0
U_global(nnodes_m*ncoord+1:end,1) = cp0;
    
elseif problem == 5
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\convergence_5.txt', 'w');
ncoord = 2;
[nodeIDs,elementIDs,nnodes,elementNodes,left_edge_nodes,bottom_edge_nodes,...
                           right_edge_nodes,top_edge_nodes,nodeCoords,~] = inp_coordinates(ncoord,problem);
nnode_ele_m = 8;
nnode_ele_c = 8;
nnodes_m = nnodes;
nnodes_c = nnodes;
numElements = length(elementIDs);
j_flux = 0;
m_d = 0.05;
cp_solvent =0;
dt = 0.01;
tend = 20;
%nt = (tend-0)/dt+1;
nt = 2001;
Lambda = linspace(1,1.35,2001);
x_coords = nodeCoords(:,1);
y_coords = nodeCoords(:,2);
stretch = ones(numElements,nt);
sigma_xx = zeros(numElements,4,nt);
sigma_yy = zeros(numElements,4,nt);
sigma_xy = zeros(numElements,4,nt);

right_top = find(abs(x_coords - 1) < 1e-6 & y_coords >= 1);
right_bottom = find(abs(x_coords - 1) < 1e-6 & y_coords < 1);



% Initialize global degree of freedom
U_global = zeros(nnodes_m*ncoord+nnodes_c,nt);  
% Initial chemical potential [μ]_0
U_global(nnodes_m*ncoord+1:end,1) = cp0;
%U_global(nnodes_m*ncoord+right_top,1) = -0.2;
%U_global(nnodes_m*ncoord+top_edge_nodes,1) = -0.2;
%U_global(nnodes_m*ncoord+right_bottom,1) = 0;
%U_global(nnodes_m*ncoord+bottom_edge_nodes,1) = 0;

end

% Material parameters
[G0,kB,NU,nu,X,md,mc,m] = material_parameters(m_d,ncoord);

time = zeros(nt,1);
max_iterations = 99;
%Time increment loop
for t =2:nt
    time(t) = time(t-1)+dt;
    if time(t)<=1
        j_cap = j_flux;
    else
        j_cap = 0;
    end
   % lambda = Lambda(t);
   % U_global(ncoord*(1:8)-1,t) = (lambda-1)*x_coords;
   % U_global(ncoord*(1:8),t) = (lambda-1)*y_coords;
   % U_global(nnodes*ncoord+1:end,t) = prescribed_cp(t);
    
    % Initialize current degree of freedom at t_n+1,[w;μ]_n+1
    U_nr = zeros(nnodes_m*ncoord+nnodes_c,max_iterations);
     % Initialize incremental degree of freedom [Δw;Δμ]
    du = zeros(nnodes_m*ncoord+nnodes_c,max_iterations);
    
    % Take converged solution from previous time step t_n as initial guess
    if t > 2
       U_nr(:, 1) = U_global(:,t-1)+(U_global(:,t-1)-U_global(:,t-2));
    else
        U_nr(:, 1) = U_global(:,t-1);
    end
 
% Newton Raphson loop
for k =1:max_iterations  
    % Initialize global stiffness matrices and global force vectors
    K_global_mat = zeros(nnodes_m*ncoord);                                 % Global material stiffness matrix
    K_global_geo= zeros(nnodes_m*ncoord);                                  % Global geometric stiffness matrix
    P_global = zeros(nnodes_m*ncoord,nnodes_c);                            % Global coupling stiffness matrix
    H_global = zeros(nnodes_c);                                            % Global chemical stiffness matrix
    H_edge_global = zeros(nnodes_c);                                       % Global chemical convection stiffness matrix
    F_int_m_global = zeros(nnodes_m*ncoord,1);                             % Global internal mechanical force vector
    F_ext_m_global = zeros(nnodes_m*ncoord,1);                             % Global external mechanical force vector
    F_int_c_global = zeros(nnodes_c,1);                                    % Global internal chemical force vector
    F_int_c_edge_global = zeros(nnodes_c,1);                               % Global internal chemical convection force vector
    F_ext_c_global = zeros(nnodes_c,1);                                    % Global external chemical force vector
% Element loop
for element_id = 1:numElements
% Element nodes and coordinates
element_nodes = elementNodes(element_id, :);
element_coordinates = nodeCoords(element_nodes, :);

% Global displacement degree of freedom for element
if problem == 3
     global_dofs = [ncoord*element_nodes-2; ncoord*element_nodes-1;ncoord*element_nodes];
else
     global_dofs = [ncoord*element_nodes-1; ncoord*element_nodes];
end
% Element nodal displacements [w_n+1]
element_displacement = U_nr(global_dofs(:), k);
 
% Element nodal increment displacements [w_n+1 - w_n]
ele_incr = U_nr(global_dofs(:), k) - U_global(global_dofs(:), t-1);

% Element nodal chemical potential [μ]_n+1
element_cp = U_nr(nnodes_m*ncoord+element_nodes(1:nnode_ele_c), k);

% Get local node indices of boundary elements on which flux or solvent is applied
if problem == 5
%corner_local_indices = [1, 2, 3, 4]; % Only these carry μ (Q8μ4)
edge_definitions = {[1,5, 2], [2,6, 3], [3,7, 4], [4,8, 1]}; % CCW edges

boundary_local_nodes = [];

for iedge = 1:4
    edge_nodes = edge_definitions{iedge};
    global_nodes = element_nodes(edge_nodes); % Get global IDs of edge nodes 
    % Check if edge is on any boundary
if all(ismember(global_nodes, bottom_edge_nodes))
        boundary_type = 'bottom';
        cp_solvent = 0;
elseif all(ismember(global_nodes, top_edge_nodes))
        boundary_type = 'top';
        cp_solvent = -0.2;
elseif all(ismember(global_nodes, right_edge_nodes))
        boundary_type = 'right';
        % Handle vertical gradient (μ_ext = -0.2 if y > 1.0, else 0)
elseif all(ismember(global_nodes, left_edge_nodes))
        boundary_type = 'left';
        continue; % No flux (natural BC)
else
        continue; % Not a boundary edge
end
boundary_local_nodes = [boundary_local_nodes, edge_nodes];
end
%boundary_local_nodes = unique(boundary_local_nodes);
else
[~, boundary_local_nodes] = ismember(flux_nodes, element_nodes);
boundary_local_nodes(boundary_local_nodes == 0) = [];
end

boundary_local_nodes;
    
    % Call element routine
    [K_ele_mat,K_ele_geo,J,F,K_stress,P_ele,H_ele,H_ele_edge,F_int_ele_m,F_int_ele_c,F_int_ele_c_edge,F_ext_ele_c,cauchy_stress,D,j] = element_routine(ncoord,nnode_ele_m,nnode_ele_c,element_coordinates,element_displacement,...
                                          ele_incr,element_cp,boundary_local_nodes,G0,kB,NU,nu,X,md,mc,m,dt,problem,element_id,j_cap,cp_solvent,TestCase,patchTest,deformControl);
    
    stretch(element_id,t) = mean(J);                   % Stretch of element for t time step
    sigma_xx(element_id,:,t) = cauchy_stress(1,1,:);   % Stress_xx for element_id at integration points for t time step
    sigma_yy(element_id,:,t) = cauchy_stress(2,2,:);   % Stress_yy for element_id at integration points for t time step
    sigma_xy(element_id,:,t) = cauchy_stress(1,2,:);   % Stress_xy for element_id at integration points for t time step
    fluid_flux(element_id,:,:,t) = j;                  % Fluid flux for element_id at integration points for t time step
    deform_drag(:,:,element_id,t) = F;                 % Deformation gradient for element_id at integration points for t time step
   
    % Assemble mechanical global stiffness matrices K_mat, K_geo
    K_global_mat(global_dofs(:),global_dofs(:)) = K_global_mat(global_dofs(:),global_dofs(:))+K_ele_mat;
    K_global_geo(global_dofs(:),global_dofs(:)) = K_global_geo(global_dofs(:),global_dofs(:))+ K_ele_geo;
    
    % Assemble global coupled stiffness matrix P
    P_global(global_dofs(:),element_nodes(1:nnode_ele_c)) = P_global(global_dofs(:),element_nodes(1:nnode_ele_c)) + P_ele;
    
    % Assemble chemical stiffness matrix H
    H_global(element_nodes(1:nnode_ele_c),element_nodes(1:nnode_ele_c)) = ...
    H_global(element_nodes(1:nnode_ele_c),element_nodes(1:nnode_ele_c)) + H_ele;
     
    % Assemble boundary chemical stiffness H_edge/surface
    H_edge_global(element_nodes(1:nnode_ele_c),element_nodes(1:nnode_ele_c)) = ...
    H_edge_global(element_nodes(1:nnode_ele_c),element_nodes(1:nnode_ele_c)) + H_ele_edge; 
                                                        
    % Assemble internal chemical forces F_int_c
    F_int_c_global(element_nodes(1:nnode_ele_c)) = F_int_c_global(element_nodes(1:nnode_ele_c)) + F_int_ele_c;
    
    % Assemble internal boundary chemical forces F_int_c_edge/surface                                               
    F_int_c_edge_global(element_nodes(1:nnode_ele_c)) = ...
    F_int_c_edge_global(element_nodes(1:nnode_ele_c)) + F_int_ele_c_edge;
                                                                  
    % Assemble external chemical forces F_ext_c
    F_ext_c_global(element_nodes(1:nnode_ele_c)) = F_ext_c_global(element_nodes(1:nnode_ele_c)) + F_ext_ele_c;                                                          
  
    % Assemble internal mechanical forces F_int_m
    F_int_m_global(global_dofs(:)) = F_int_m_global(global_dofs(:)) + F_int_ele_m;
    
    %Since there are no body forces and no tractions F_ext_m_global = 0
    
end % End of element loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PATCH TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(patchTest, 'Yes')
   f = 1000;
   k_global = K_global_mat;
   f_int = F_int_m_global;
   %boundary conditions
   %both x and y displacements fixed at node 1
   k_global(1,:) = 0; 
   k_global(:,1) = 0;
   k_global(2,:) = 0; 
   k_global(:,2) = 0;
   k_global(1,1) = 1;
   k_global(2,2) = 1;
   %x displacements at node 4 are fixed
   k_global(7,:) = 0;
   k_global(:,7) = 0;
   k_global(7,7) = 1;
   % applied forces
   f_int(5) = f;
   f_int(11) = 2*f;
   f_int(13) = -f;
   f_int(17) = f;

  tu = k_global\f_int; % Solve for diaplacements
  
  strain_xx = zeros();
  strain_yy = zeros();
  strain_xy = zeros();
  strain_yx = zeros();
  stress_xx = zeros();
  stress_yy = zeros();
  stress_xy = zeros();
  stress_yx = zeros();
  figure;
  hold on;
  axis equal;
  % Post processing, compute stress, strain, tests and plot
  for element = 1:numElements
    element_nodes = elementNodes(element, :);                                     % Element nodes
    element_coordinates = nodeCoords(element_nodes, :);                           % Coordinates of the element
    expected_area = polyarea(element_coordinates(:,1),element_coordinates(:,2));  % Area of the irregular element
    global_dofs = [ncoord*element_nodes-1; ncoord*element_nodes];                 % Global nodal displacement dof
    element_displacement = tu(global_dofs(:));
    coords = [element_coordinates; element_coordinates(1,:)];                     % close the polygon to plot
    %plot(coords(:,1), coords(:,2), 'k-');
    fill(coords([1:end, 1], 1), coords([1:end, 1], 2), 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'blue');
    % Element number at centroid
    centroid = mean(nodeCoords(element_nodes, :), 1);
    text(centroid(1), centroid(2), num2str(element), 'Color', 'b', 'FontSize', 12);
    actual_area = 0;
    for ip = 1:4
        [xi,eta,zeta,w] = integration_points(ip,2);
        [N_u,N_cp] = shape_functions(xi,eta,zeta,1);                                     
        [dN_u,dN_cp] = shape_function_derivatives(xi,eta,zeta,1); 
        B_u = zeros(3,8);
        % Jacobian matrix
        Jacobian = dN_u*element_coordinates; detjacobian = det(Jacobian);
        invjacobian = inv(Jacobian);
        actual_area = actual_area + detjacobian;
        B = invjacobian*dN_u;
        for a = 1:4
        B_u(1,2*a-1) = B(1,a); % B matrix
        B_u(2,2*a-1) = 0;
        B_u(3,2*a-1) = B(2,a);
        B_u(1,2*a)   = 0;
        B_u(2,2*a)   = B(2,a);
        B_u(3,2*a)   = B(1,a);
        end
    % Compute strain
    strain = B_u*element_displacement;
    % Compute stress
    stress = D*strain;
    strain_xx(element,ip) = strain(1);
    strain_yy(element,ip) = strain(2);
    strain_xy(element,ip) = strain(3)/2;
    strain_yx(element,ip) = strain(3)/2;
    stress_xx(element,ip) = stress(1);
    stress_yy(element,ip) = stress(2);
    stress_xy(element,ip) = stress(3);
    stress_yx(element,ip) = stress(3);
    end
   % Area of element calculated from determinant of jacobian test
   if (abs(expected_area)- actual_area<1e-10)
      cprintf('*Comments', 'Element id: %d\n', element);
      fprintf('Actual area : %f \nArea from jacobian : %f \n', actual_area, expected_area);
      cprintf('blue','Area test passed\n');
   end
  end
% Test if stress or strain at integration points is same
% Change strain_xx to strain_yy, stress_xx to repective parameter
all_equal = all(abs(stress_xx(:))- stress_xx(1)<1e-10);
if all_equal
fprintf('Stress in xx direction : \n'); disp(stress_xx);
fprintf('Strain in xx direction : \n'); disp(strain_xx);
cprintf('blue','Strain and Stress is equal at all integration points\n');
cprintf('blue','Patch test passed\n');
end 
  
  % Plot node numbers
  numNodes = size(nodeCoords, 1);
  for n = 1:numNodes
  x = nodeCoords(n,1);
  y = nodeCoords(n,2);
  plot(x, y, 'ro','MarkerSize', 6, 'MarkerFaceColor', 'r');  % node marker
  text(x - 0.09, y - 0.07, num2str(n), 'Color', 'k', 'FontSize', 12)
  
  f_v = 1;   
  arrowScale = 0.2;
% Example force vectors: [x, y] per node
% Adjust indices (e.g., nodeCoords(n, :)) according to your model
forceVectors = [ f_v,  0;   2*f_v,  0;   f_v,  0;   -f_v, 0];
forceNodes = [3, 6, 9, 7]; 
% Plot arrows
for i = 1:length(forceNodes)
    node = forceNodes(i);
    x = nodeCoords(node, 1);
    y = nodeCoords(node, 2);
    u = forceVectors(i, 1);  % X-component of force
    v = forceVectors(i, 2);  % Y-component of force
    quiver(x, y, arrowScale*u, arrowScale*v, 0, ...
        'LineWidth', 2, 'MaxHeadSize', 1, 'Color', 'k');
    labels = {'F', '2F', 'F', 'F'};
    text(x + 0.1, y + 0.09, labels{i}, 'FontSize', 12, 'FontWeight', 'bold');

end
  
  end
  %title('Finite Element Patch');
  xlabel('X'); ylabel('Y'); xlim([-0.5 3]);ylim([-1 2.5])
 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form coupled stiffness matrix and residual
H_global = H_global - H_edge_global ;
F_int_c_global = F_int_c_global - F_int_c_edge_global;
K_global =  K_global_mat + K_global_geo;
K_C_global = [K_global,P_global;P_global.',H_global];
F_int = [F_int_m_global; F_int_c_global];
F_ext = [F_ext_m_global; F_ext_c_global];
R_global = F_int - F_ext ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
if problem == 1
%Y displacements are constrained for bottom edge
y_fixed = 2*bottom_edge_nodes.';
K_C_global(y_fixed,:) = 0;
K_C_global(:,y_fixed) = 0;
K_C_global(sub2ind(size(K_C_global),y_fixed,y_fixed)) = 1;
R_global(y_fixed) = 0;
% X displacements are constrained for all nodes to solve for 1D problem
x_fixed = 2*(1:nnodes)-1;
K_C_global(x_fixed,:) = 0;
K_C_global(:,x_fixed) = 0;
K_C_global(sub2ind(size(K_C_global),x_fixed,x_fixed)) = 1;
R_global(x_fixed) = 0;
% Y displacements are prescribed and constrained for top edge
y_top_dof = ncoord*flux_nodes;
K_C_global(y_top_dof,:) = 0;
K_C_global(:,y_top_dof) = 0;
K_C_global(sub2ind(size(K_C_global),y_top_dof,y_top_dof)) = 1;
R_global(y_top_dof) = 0;
% Chemical potential at top is made equal to cp_solvent and fixed
cp_fixed = nnodes*ncoord+flux_nodes;
K_C_global(cp_fixed,:) = 0;
K_C_global(:,cp_fixed) = 0;
K_C_global(sub2ind(size(K_C_global),cp_fixed,cp_fixed)) = 1;
R_global(cp_fixed) = 0;

elseif problem == 2 
% Different BC for deformation control
if strcmpi(deformControl,'Yes')
% X displacements are constrained for left edge
bdof = 2*[1,4].'-1;
K_C_global(bdof,:) = 0;
K_C_global(:,bdof) = 0;
K_C_global(sub2ind(size(K_C_global),bdof,bdof)) = 1;
R_global(bdof) = 0;
% Y displacement is fixed node 2
bnode = ncoord*2;
K_C_global(bnode,:) = 0;
K_C_global(:,bnode) = 0;
K_C_global(bnode,bnode) = 1;
R_global(bnode) = 0;
else
% X displacements are constrained for AD edge
bdof = 2*left_edge_nodes.'-1;
K_C_global(bdof,:) = 0;
K_C_global(:,bdof) = 0;
K_C_global(sub2ind(size(K_C_global),bdof,bdof)) = 1;
R_global(bdof) = 0;
% Y displacement is constrained at B (node 6)
bnode = ncoord*6;
K_C_global(bnode,:) = 0;
K_C_global(:,bnode) = 0;
K_C_global(bnode,bnode) = 1;
R_global(bnode) = 0;
end

elseif problem == 3
% Y displacements are constrained on the plane xoz
xoz_fixed_dof = ncoord*xoz_fixed- 1;
K_C_global(xoz_fixed_dof,:) = 0;
K_C_global(:,xoz_fixed_dof) = 0;
K_C_global(sub2ind(size(K_C_global),xoz_fixed_dof,xoz_fixed_dof)) = 1;
R_global(xoz_fixed_dof) = 0;
% X displacements are constrained on the plane yoz 
yoz_fixed_dof = ncoord*yoz_fixed - 2;
K_C_global(yoz_fixed_dof,:) = 0;
K_C_global(:,yoz_fixed_dof) = 0;
K_C_global(sub2ind(size(K_C_global),yoz_fixed_dof,yoz_fixed_dof)) = 1;
R_global(yoz_fixed_dof) = 0;
% Z displacements are constrained for point C
z_fixed = ncoord*z_fixed_node;
K_C_global(z_fixed,:) = 0;
K_C_global(:,z_fixed) = 0;
K_C_global(sub2ind(size(K_C_global),z_fixed,z_fixed)) = 1;
R_global(z_fixed) = 0;

elseif problem == 4
%Symmetrical boundary conditions on left and bottom edges
% X displacements are constrained for AD edge
left_dof = 2*left_edge_nodes.'-1;
K_C_global(left_dof,:) = 0;
K_C_global(:,left_dof) = 0;
K_C_global(sub2ind(size(K_C_global),left_dof,left_dof)) = 1;
R_global(left_dof) = 0; 
% Y displacements are constrained for AB edge
bottom_dof = 2*bottom_edge_nodes.';
K_C_global(bottom_dof,:) = 0;
K_C_global(:,bottom_dof) = 0;
K_C_global(sub2ind(size(K_C_global),bottom_dof,bottom_dof)) = 1;
R_global(bottom_dof) = 0;

elseif problem == 5
%Symmetrical boundary conditions on left edge
% X displacements are constrained for AD edge
left_dof = 2*left_edge_nodes.'-1;
K_C_global(left_dof,:) = 0;
K_C_global(:,left_dof) = 0;
K_C_global(sub2ind(size(K_C_global),left_dof,left_dof)) = 1;
R_global(left_dof) = 0;

% Y displacement is constrained at point E 
fixed_nodes = 36;
fixed_dof = 2*fixed_nodes;
K_C_global(fixed_dof,:) = 0;
K_C_global(:,fixed_dof) = 0;
K_C_global(sub2ind(size(K_C_global),fixed_dof,fixed_dof)) = 1;
R_global(fixed_dof) = 0;

%{
right_top_fix = 2*nnodes+right_top;
K_C_global(right_top_fix,:) = 0;
K_C_global(:,right_top_fix) = 0;
K_C_global(sub2ind(size(K_C_global),right_top_fix,right_top_fix)) = 1;
R_global(right_top_fix) = 0;

right_bottom_fix = 2*nnodes+right_bottom;
K_C_global(right_bottom_fix,:) = 0;
K_C_global(:,right_bottom_fix) = 0;
K_C_global(sub2ind(size(K_C_global),right_bottom_fix,right_bottom_fix)) = 1;
R_global(right_bottom_fix) = 0;

top_fix = 2*nnodes+top_edge_nodes;
K_C_global(top_fix,:) = 0;
K_C_global(:,top_fix) = 0;
K_C_global(sub2ind(size(K_C_global),top_fix,top_fix)) = 1;
R_global(top_fix) = 0;

bottom_fix = 2*nnodes+bottom_edge_nodes;
K_C_global(bottom_fix,:) = 0;
K_C_global(:,bottom_fix) = 0;
K_C_global(sub2ind(size(K_C_global),bottom_fix,bottom_fix)) = 1;
R_global(bottom_fix) = 0;
%}

end

res(k,t)= norm(R_global);

K_inv = inv(K_C_global);

% Solve for incremental degere of freedom Δu = [Δw;Δμ]
du(:,k) = -K_C_global\R_global;

% Update Newton Raphson iteration [w;μ]_k+1 = [w;μ]_k + [Δw;Δμ]
if problem == 5
U_nr(:,k+1) = U_nr(:,k)+ 1*du(:,k);
else
U_nr(:,k+1) = U_nr(:,k)+  du(:,k);
end

err_du(k,t) =  norm(du(1:ncoord*nnodes,k));
err_dcp(k,t) =  norm(du(ncoord*nnodes+1:end,k));
res_u(k,t) = norm(R_global(1:ncoord*nnodes));
res_cp(k,t) = norm(R_global(ncoord*nnodes+1:end));
error_u(k,t) = norm(du(1:ncoord*nnodes,k))/norm(U_nr(1:ncoord*nnodes,k+1));
error_cp(k,t) = norm(du(ncoord*nnodes+1:end,k))/norm(U_nr(ncoord*nnodes+1:end,k+1));

error(k,t) = norm(du(:,k))/norm(U_nr(:,k+1));
if error(k,t) < 1e-6
    fprintf('Convergence at iteration %d\n for time step %d\n', k,t);
    % Write to a file for postprocessing
    fprintf(fileID, 'Convergence achieved at %d iteration for time step %d\n', k,t);
    break;
end

end %End of Newton Raphson loop

%Update global degree of freedom
U_global(:,t) =  U_nr(:,k+1);

end % End of time loop
fclose(fileID);
% Nodal chemical potential
cp_global = U_global(nnodes*ncoord+1:end,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store displacement vector, chemical potential, stress for plotting
if problem == 1 && strcmpi(patchTest,'No')

% Evaluate chemical potential from stretch for the element to plot
mu = (G0*nu.*stretch - G0*nu*1./stretch+...
                        kB*NU*(log(1-0.999./stretch)+1./stretch + X./stretch.^2))/(kB*NU);
                    
% Evaluate stress from stretch and chemical potential
sigma_x = -kB*NU/nu.*mu + kB*NU/nu*(log(1-0.999./stretch)+1./stretch + X./stretch.^2);

% Store stretch
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\stretch.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d elements, %d time steps\n', numElements, nt);
fclose(fileID); 
dlmwrite(filename, stretch, '-append', 'delimiter', '\t', 'precision', '%e');

% Store chemical potential
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\chemical_potential.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d elements, %d time steps\n', numElements, nt);
fclose(fileID); 
dlmwrite(filename, mu, '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\stress.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d elements, %d time steps\n', numElements, nt);
fclose(fileID); 
dlmwrite(filename, sigma_x, '-append', 'delimiter', '\t', 'precision', '%e');

% Store error for all time steps
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\error.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d iterations, %d time steps\n', max_iterations, nt);
fclose(fileID); 
dlmwrite(filename, error, '-append', 'delimiter', '\t', 'precision', '%e');

elseif problem == 2 && strcmpi(deformControl,'No')

% Store displacement vector
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\displacement_vector.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Displacement vector: %d nodes, %d DOFs, %d time steps\n', nnodes, ncoord*nnodes, nt);
fclose(fileID); 
dlmwrite(filename, U_global(1:nnodes*ncoord,:), '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress in xx direction
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\stress_xx.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stress at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, sigma_xx, '-append', 'delimiter', '\t', 'precision', '%e');

% Store chemical potential
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\chemical_potential.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Chemical potential at nodes : %d nodes, %d time steps\n', ...
        numElements,nt);
fclose(fileID); 
dlmwrite(filename, cp_global, '-append', 'delimiter', '\t', 'precision', '%e');

% Store error
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\error.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d iterations, %d time steps\n', max_iterations, nt);
fclose(fileID); 
dlmwrite(filename, error, '-append', 'delimiter', '\t', 'precision', '%e');

elseif problem == 3
    
% Store displacements
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\displacement_vector.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Displacement vector: %d nodes, %d DOFs, %d time steps\n', nnodes, ncoord*nnodes, nt);
fclose(fileID); 
dlmwrite(filename, U_global(1:nnodes*ncoord,:), '-append', 'delimiter', '\t', 'precision', '%e');

% Store error
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\error3.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d iterations, %d time steps\n', max_iterations, nt);
fclose(fileID); 
dlmwrite(filename, error, '-append', 'delimiter', '\t', 'precision', '%e');


elseif problem == 4
    
% Store displacements
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\displacement_vector.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Displacement vector: %d nodes, %d DOFs, %d time steps\n', nnodes, ncoord*nnodes, nt);
fclose(fileID); 
dlmwrite(filename, U_global(1:nnodes*ncoord,:), '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress in xx direction
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\stress_xx.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stress at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, sigma_xx, '-append', 'delimiter', '\t', 'precision', '%e');

% Store nodal chemical potential
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\chemical_potential.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Chemical potential at nodes : %d nodes, %d time steps\n', ...
        numElements,nt);
fclose(fileID); 
dlmwrite(filename, cp_global, '-append', 'delimiter', '\t', 'precision', '%e');

% Store error
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\error4.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d iterations, %d time steps\n', max_iterations, nt);
fclose(fileID); 
dlmwrite(filename, error, '-append', 'delimiter', '\t', 'precision', '%e');

elseif problem == 5
 % Store displacements
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\displacement_vector.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Displacement vector: %d nodes, %d DOFs, %d time steps\n', nnodes, ncoord*nnodes, nt);
fclose(fileID); 
dlmwrite(filename, U_global(1:nnodes*ncoord,:), '-append', 'delimiter', '\t', 'precision', '%e');

% Store nodal chemical potential
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\chemical_potential.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Chemical potential at nodes : %d nodes, %d time steps\n', ...
        numElements,nt);
fclose(fileID); 
dlmwrite(filename, cp_global, '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress in xx direction
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\stress_xx.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stress at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, sigma_xx, '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress in yy direction
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\stress_yy.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stress at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, sigma_yy, '-append', 'delimiter', '\t', 'precision', '%e');

% Store stress in xy direction
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\stress_xy.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stress at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, sigma_xy, '-append', 'delimiter', '\t', 'precision', '%e');

% Store fluid flux
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\fluid_flux.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Flux at Gauss points: %d elements, %d points, %d time steps\n', ...
        numElements, 4, nt);
fclose(fileID); 
dlmwrite(filename, fluid_flux, '-append', 'delimiter', '\t', 'precision', '%e');

% Store error
filename = "C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\error5.dat";
fileID = fopen(filename, 'w');
fprintf(fileID, '# Stretch: %d iterations, %d time steps\n', max_iterations, nt);
fclose(fileID); 
dlmwrite(filename, error, '-append', 'delimiter', '\t', 'precision', '%e');
end

% Call plots function
plots_function(problem);