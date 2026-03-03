%% Calculates stiffness matrices and force vectors of an element
%__________________________________________________________________________
% INPUT: Model and material parameters
%__________________________________________________________________________
% OUTPUT: Element stiffness matrices, force vectors, cauchy stress
%__________________________________________________________________________
function [K_ele_mat,K_ele_geo,J,F,K_stress,P_ele,H_ele,H_ele_edge,F_int_ele_m,F_int_ele_c,F_int_ele_c_edge,F_ext_ele_c,cauchy_stress,D,j] = element_routine(ncoord,nnode_ele_m,nnode_ele_c,element_coordinates,element_displacement,...
                                          ele_incr,element_cp,boundary_local_nodes,G0,kB,NU,nu,X,~,mc,m,dt,problem,element_id,j_cap,cp_solvent,TestCase,patchTest,deformControl)
   if ncoord == 2
      nip = 4; 
   elseif ncoord ==3
      nip = 8;
   end
   % Initialize element stiffness matrices and force vectors
   K_ele_mat = zeros(ncoord*nnode_ele_m);           % Element material stiffness matrix
   K_ele_geo = zeros(ncoord*nnode_ele_m);           % Element geometric stiffness matrix
   P_ele = zeros(ncoord*nnode_ele_m,nnode_ele_c);   % Element coupling stiffness matrix
   H_ele = zeros(nnode_ele_c);                      % Element chemical stiffness matrix
   H_ele_edge = zeros(nnode_ele_c);                 % Element chemical stiffness matrix,convection term               
   F_int_ele_m = zeros(ncoord*nnode_ele_m,1);       % Internal mechanical force matrix
   F_int_ele_c = zeros(nnode_ele_c,1);              % Internal chemical force matrix
   F_int_ele_c_edge = zeros(nnode_ele_c,1);         % Internal chemical force vector,convection term
   F_ext_ele_c = zeros(nnode_ele_c,1);              % External chemical force matrix
   J = zeros(nip,1);                                % Determinant of deformation gradient
   cauchy_stress = zeros(ncoord,ncoord,nip);        % Cauchy stress at integration points
   K_stress = zeros(ncoord,ncoord,nip);             % Kirchoff stress at integration points
   area = 0;                                        % Initialize area to accumulate determinant of jacobian
   for ip = 1:nip % Integration loop
       [xi,eta,zeta,w] = integration_points(ip,ncoord);                                          % Integration points and weights
       [~,N_cp] = shape_functions(xi,eta,zeta,problem);                                          % Shape functions
       [dN_u,dN_cp] = shape_function_derivatives(xi,eta,zeta,problem);                           % Derivatives of shape functions
       [~,invjacobian,detjacobian,invjacob_current,~,invjacob_current_c]= ...                    % Jacobian,determinant of jacobian,inverse jacobian
           jacobian(nnode_ele_c,dN_u,dN_cp,element_coordinates,element_displacement);   
       [G,G_grad,B_u,B_cp,B_p,B_f] = B_matrices(invjacobian,invjacob_current,...                 % B matrices
           invjacob_current_c,dN_u,dN_cp,ncoord,nnode_ele_m,nnode_ele_c,problem);     
       % For problem == 1, prescribe displacements for the boundary element
       % and for TestCase = 'Yes' no need to prescribe displacements
       if problem == 1&& element_id == 20 
           if strcmpi(TestCase, 'No') &&  strcmpi(patchTest, 'No')
           element_displacement = [0;0.473099994024000;0;0.473099994024000;...
                                   0;0.498000000000000;0;0.498000000000000];
           end
       end 
       % Displacement gradient [dNi/dx,dNi/dy]*[ux,uy]
       % Deformatiom gradient F = I + Grad(u)    
       if ncoord == 3
       displacement_gradient =  G_grad*element_displacement;   
       F = eye(ncoord) + [displacement_gradient(1), displacement_gradient(2), displacement_gradient(3);
                          displacement_gradient(4), displacement_gradient(5), displacement_gradient(6);
                          displacement_gradient(7), displacement_gradient(8), displacement_gradient(9)];
       else
       displacement_gradient =  G_grad*element_displacement;
       F = eye(ncoord) + [displacement_gradient(1), displacement_gradient(2);...
                          displacement_gradient(3), displacement_gradient(4)];                
       end
       % Determinant of deformation gradient
       J(ip) = det(F);  
       % Chemical potential of an element
       cpe = N_cp*element_cp;   
       % Material routine for material tangent modulus, Cauchy stress and Kirchoff stress at integration points
       [D, C_stress, K_stress(:,:,ip),S]= material_routine(F,J(ip),cpe,G0,kB,NU,nu,X,ncoord);
       % Store cauchy stress and flux for plotting
       cauchy_stress(:,:,ip) = C_stress;
       j(:,:,ip) = - kB*NU*kB*NU*m*B_cp*element_cp; % j = -mgradμ
       
       if problem ==1 && strcmpi(patchTest, 'Yes')
         % Define material properties for patch test
         E = 1e7;     % Young's modulus(Pa)
         nu = 0.3;    % Poisson's ratio
       % Plane strain D matrix
        D = (E / ((1 + nu) * (1 - 2 * nu))) * ...
                  [1 - nu,   nu,         0;
                  nu, 1 - nu,     0;
                  0, 0, (1 - 2 * nu) / 2];
       end
       % Voigt notation
       if ncoord == 2
       K_stress_v = [K_stress(1,1,ip);K_stress(2,2,ip); K_stress(1,2,ip)];          % τ_v(3x1) = [τ_xx ; τ_xy ; τ_yx], B_u(3x16)
       K_stress_hat = blkdiag(K_stress(:,:,ip),K_stress(:,:,ip));                   % τ_hat(4x4) = diag(τ,τ) with  τ = [τ_xx,τ_xy ; τ_yx,τ_yy]
       elseif ncoord == 3
       K_stress_v = [K_stress(1,1,ip); K_stress(2,2,ip); K_stress(3,3,ip); ...      % τ_v(6x1) = [τ_xx ; τ_yy ; τ_zz ; τ_xy ; τ_xz ; τ_yz]
                     K_stress(1,2,ip); K_stress(1,3,ip); K_stress(2,3,ip)];
       K_stress_hat = blkdiag(K_stress(:,:,ip),K_stress(:,:,ip),K_stress(:,:,ip));  % τ_hat(9x9) = diag(τ,τ,τ)                                                                                                              %   τ_yx, τ_yy , τ_yz;             
       end                                                                                                            
       
       K_ele_mat = K_ele_mat + (w*detjacobian*B_u.'*D*B_u);                                 % Material stiffness matrix eq.31(1)      
       K_ele_geo = K_ele_geo + (w*detjacobian*G.'*K_stress_hat*G);                          % Geometric stiffness matrix eq.31(2)
       P_ele = P_ele + w*detjacobian*(- J(ip)/nu*kB*NU*B_p.'*N_cp);                         % Coupling stiffness matrx eq.32(1)
       H_ele = H_ele + w*detjacobian*(- J(ip)*kB*NU*kB*NU*dt*B_cp.'*m*B_cp);                % Chemical stiffness matrix eq.32(2)
       F_int_ele_m = F_int_ele_m + w*detjacobian*(B_u.'*K_stress_v);                        %Internal mechanical force vector eq.27(1) 
       F_int_ele_c = F_int_ele_c + w*detjacobian*(- J(ip)/nu*kB*NU*N_cp.'*B_p*ele_incr...   %Internal chemical force vector eq.27(2)
                                        - J(ip)*kB*NU*kB*NU*dt*B_cp.'*m*B_cp*element_cp);  
       % Calculate area form determinant of jacobian for test case                           
       area =  area + detjacobian; 
%_________________________________________________________________________
 %% Test tangent modulus
 %% Aim: To compute finite difference approximation tangent modulus and compare to FEM tangen modulus D
 %verify_tangent_modulus(F,cpe,ncoord)  % Uncomment for tangent test
   end % End of integration loop
   
%% %%%%%%%%%%%%%%%%%%%%%%% Boundary terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The convection/boundary terms F_int_ele_c_edge, H_ele_edge and external
% chemical force vector F_ext_ele_c are to be evaluated on the boundary,
% 1D integrtaion over the edge of an element for 2D elements and surface 2D
% integration for 3D elements on the nodes where fluid flux or external solvent is applied 
 if problem == 1
     l_e = 0;
     boundary_local = boundary_local_nodes;
      if  ~isempty(boundary_local) && length(boundary_local)>1
          boundary_coordinates = element_coordinates(boundary_local,:);
          if boundary_coordinates(:,2) == 1
          %Defining local coordinate values for a 2 point integration
          xip = [-1/sqrt(3),1/sqrt(3)]; etap = [1,1]; zeta = 0;% Gauss points
          wp = [1,1];  % Weights
          n_local = 1; % n_local = 1 is [N,xi] for the nodes on top surface 
          end
          %Integrationg on the edge
          for i = 1:2
          xi = xip(i);
          eta = etap(i);
          [~,N_cp_edge] = shape_functions(xi,eta,zeta,problem);
          [dN_edge] = shape_function_derivatives(xi,eta,zeta,problem);
          jacob = norm(dN_edge(n_local,:)*(element_coordinates));
          l_e = l_e+jacob;
          F_ext_ele_c = F_ext_ele_c - kB*NU*kB*NU*dt*mc*wp(i)*jacob*(N_cp_edge.'*cp_solvent);
          end
          T = l_e/6*[2,1;1,2]; % Taken from reference paper
          H_ele_edge(boundary_local,boundary_local) =  kB*NU*kB*NU*dt*mc*T;
          F_int_ele_c_edge(boundary_local) =  kB*NU*kB*NU*dt*mc*T*element_cp(boundary_local);
      end
       
  elseif problem == 2
     boundary_local = boundary_local_nodes;
     l_e = 0;
      if ~isempty(boundary_local) && length(boundary_local)>1
          boundary_coordinates = element_coordinates(boundary_local,:);
          if boundary_coordinates(:,2) == 1
          %Defining local coordinate values for a 2 point integration
          xip = [-1/sqrt(3),1/sqrt(3)]; etap = [1,1]; zeta = 0;% Gauss points
          wp = [1,1]; % Weights
          n_local = 1;
          elseif boundary_coordinates(:,1)==1
          etap = [-1/sqrt(3),1/sqrt(3)]; xip = [1,1]; zeta = 0;% Gauss points
          wp = [1,1]; % Weights
          n_local = 2;
          end
          %Integrationg on the edge
          for i = 1:2
          xi = xip(i);
          eta = etap(i);
          [~,N_cp_edge] = shape_functions(xi,eta,zeta,problem);
          [dN_edge] = shape_function_derivatives(xi,eta,zeta,problem);
          jacob = norm(dN_edge(n_local,:)*(element_coordinates));
          l_e = l_e + jacob;
          F_ext_ele_c = F_ext_ele_c+ dt*kB*NU*wp(i)*jacob*N_cp_edge.'*j_cap;
          end
         % Uncomment for deformTest = 'Yes'
         % T = l_e/30*[4,2,-1;2,16,2;-1,2,4]; % Taken from reference paper
         % H_ele_edge(boundary_local,boundary_local) =  kB*NU*kB*NU*dt*mc*T;
         % F_int_ele_c_edge(boundary_local) =  kB*NU*kB*NU*dt*mc*T*element_cp(boundary_local);
     end
    
   elseif problem == 3
      surface_area = 0;
      boundary_local = boundary_local_nodes;
     if ~isempty(boundary_local) && length(boundary_local) == 4
          boundary_coordinates = element_coordinates(boundary_local, :);
          apply_flux_bc = false;
          if all(abs(boundary_coordinates(:, 3) - 1) < 1e-6)
          x_coords = boundary_coordinates(:, 1);
          y_coords = boundary_coordinates(:, 2);
          if all(x_coords >= -1e-6 & x_coords <= 0.600000024) && all(y_coords >= -1e-6 & y_coords <= 0.600000024)
          apply_flux_bc = true;
          end
          end
          if apply_flux_bc
         % Reorder nodes counterclockwise based on x and y coordinates
          [~, sort_idx] = sortrows(boundary_coordinates(:, 1:2), [1, 2]); % Sort by x, then y
          correct_order = [1, 2, 4, 3];
          boundary_local_sorted = boundary_local(sort_idx(correct_order));
          boundary_coordinates_sorted = element_coordinates(boundary_local_sorted, :);
          xip = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3), -1/sqrt(3)];
          etap = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)];
          W = [1, 1, 1, 1];
          for i = 1:4
          xi = xip(i); eta = etap(i);
          N_cp_edge = [(1/4)*(1-xi)*(1-eta), (1/4)*(1+xi)*(1-eta), ...
                      (1/4)*(1+xi)*(1+eta), (1/4)*(1-xi)*(1+eta)];
          dN_edge = [-(1/4)*(1-eta), (1/4)*(1-eta), (1/4)*(1+eta), -(1/4)*(1+eta); ...
                     -(1/4)*(1-xi), -(1/4)*(1+xi), (1/4)*(1+xi), (1/4)*(1-xi)];
          r_xi = dN_edge(1, :) * boundary_coordinates_sorted;
          r_eta = dN_edge(2, :) * boundary_coordinates_sorted;
          jacob_vec = cross(r_xi, r_eta);
          jacob = norm(jacob_vec);
          %fprintf('Gauss point %d: jacob = %f\n', i, jacob);
          surface_area = surface_area + W(i) * jacob;
          F_ext_ele_c(boundary_local_sorted) = F_ext_ele_c(boundary_local_sorted) + kB*NU*dt* W(i)*jacob*N_cp_edge.'*j_cap;
         end
         % fprintf('Boundary face area for element %d: %f\n', element_id, area);
          end
     end 
 % For the top edge η is constant and only ξ changes where as
% for right edge η changes and ξ is constant.
  elseif problem == 4
      if ~isempty(boundary_local_nodes)
          num_edges = length(boundary_local_nodes)/2;
          for ik = 1:num_edges
          idx = (ik-1) * 2 + 1;
          l_e = 0;
          element_displacement = reshape(element_displacement, 2, 8).';
          boundary_local = boundary_local_nodes(idx:idx+1);
          boundary_coords = element_coordinates(boundary_local,:); 
          y_coords = boundary_coords(:,2);
          x_coords = boundary_coords(:,1);
          xip = [-1/sqrt(3),1/sqrt(3)]; w = 1;
          apply_bc = false;
          % solvent-exposed boundaries
          if all(abs(y_coords - 1) < 1e-6)      % Top edge
          apply_bc = true;
          elseif all(abs(x_coords - 1) < 1e-6)  % Right edge
          apply_bc = true;
          end
          if apply_bc
          % Edge length
          l_e = norm(boundary_coords(2,:) - boundary_coords(1,:)); % Edge length for test case
          L0 = 0;
          for i = 1:2
              xi = xip(i);
              N_cp = [0.5*(1 - xi), 0.5*(xi + 1)];                                  
              dN = [-0.5, 0.5];
              jacob = dN *(boundary_coords + element_displacement(boundary_local,:));
              detjacob = norm(jacob);
              L0 = L0 + w*detjacob;    
              F_ext_ele_c(boundary_local) = F_ext_ele_c(boundary_local) -...
                                              detjacob*w*dt*mc*kB*NU*kB*NU*N_cp.'*cp_solvent;   
              H_ele_edge(boundary_local, boundary_local) = H_ele_edge(boundary_local, boundary_local)+...
                                                        w*detjacob*dt*mc*kB*NU*kB*NU *(N_cp.'*N_cp);
                                                    
              F_int_ele_c_edge(boundary_local) =  F_int_ele_c_edge(boundary_local) +...
                                              w*detjacob*dt*mc*kB*NU*kB*NU *(N_cp.'*N_cp)*element_cp(boundary_local);   
          end
          % Analytical boundary term
          T = L0/6 * [2 1; 1 2];                                                       
         % H_ele_edge(boundary_local, boundary_local) = H_ele_edge(boundary_local, boundary_local)+...
         %                                              dt*mc*kB*NU*kB*NU *T;
         % F_int_ele_c_edge(boundary_local) =  F_int_ele_c_edge(boundary_local) +...
         %                                    dt*mc*kB*NU*kB*NU *T*element_cp(boundary_local);   
           end
          end
      end
% For the top and bottom edges η is constant and only ξ changes where as
% for right edge η changes and ξ is constant.
elseif problem == 5  
      element_displacement = reshape(element_displacement,2,8).';
      if ~isempty(boundary_local_nodes)
        num_edges = length(boundary_local_nodes)/3;
        for ik = 1:num_edges
        idx = (ik-1)*3 + 1;
        boundary_local = boundary_local_nodes(idx:idx+2);
        boundary_coords = element_coordinates(boundary_local,:);
        l_e = norm(boundary_coords(3,:) - boundary_coords(1,:)); % Edge length for test case
        y_coords = boundary_coords(:,2);
        x_coords = boundary_coords(:,1);
        y_center = mean(y_coords);
        xip = [-1/sqrt(3),1/sqrt(3)]; w = 1;
        apply_bc = false;  
        % Bottom edge (μ=0, allow flux)
        if all(abs(y_coords - 0) < 1e-6) %solvent_region = 'bottom';
            apply_bc = true;
            mu_ext = 0;         
        elseif all(abs(y_coords - 2) < 1e-6) % solvent_region = 'top';
            apply_bc = true;
            mu_ext = -0.2;
        elseif all(abs(x_coords - 1) < 1e-6) %solvent_region = 'right';
            apply_bc = true; 
            if y_center >1
              mu_ext = -0.2;
            else
             mu_ext = 0;
           end
        end       
        if apply_bc
        L0 =0;
        for i = 1:2
        xi = xip(i);
        N_cp = [0.5 * (xi - 1) * xi, (1 - xi) * (1 + xi), 0.5 * (xi + 1) * xi];                                  
        dN = [xi - 0.5, -2 * xi, xi + 0.5 ];
        jacob = dN *(boundary_coords + element_displacement(boundary_local,:));
        detjacob = norm(jacob);
        L0 = L0 + detjacob;
        F_ext_ele_c(boundary_local) = F_ext_ele_c(boundary_local)- w*detjacob*kB*NU*kB*NU*dt*mc*N_cp.'*mu_ext;
        F_int_ele_c_edge(boundary_local) = F_int_ele_c_edge(boundary_local) +...
                                              w*detjacob*kB*NU*kB*NU*dt*mc*(N_cp.'*N_cp)*element_cp(boundary_local);
        H_ele_edge(boundary_local, boundary_local) = H_ele_edge(boundary_local, boundary_local)+...
                                                             w*detjacob*kB*NU*kB*NU*dt*mc*(N_cp.'*N_cp);
        end
        T = L0/30 * [4, 2, -1; 2, 16, 2; -1, 2, 4];
        %F_int_ele_c_edge(boundary_local) = F_int_ele_c_edge(boundary_local) + kB*NU*kB*NU*dt*mc*T*element_cp(boundary_local);
        %H_ele_edge(boundary_local, boundary_local) = H_ele_edge(boundary_local, boundary_local)+ kB*NU*kB*NU*dt*mc*T;
        end
        end
      end
 end   % End of boundary integration terms
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runs when TestCase = 'Yes'
if strcmpi(TestCase, 'Yes')
cprintf('*Comments', 'Element id: %d\n', element_id);
testCase = matlab.unittest.TestCase.forInteractiveUse;
% Provide local coordinate values anywhere in the local domain [-1,1]
xi = -1/3; eta = -1/3 ; zeta = 1/2;
[N_u,~] = shape_functions(xi,eta,zeta,problem);                                       % Shape functions
[dN_u,~] = shape_function_derivatives(xi,eta,zeta,problem);                           % Derivatives of shape functions
%% Test case: 1          
%% Aim: To check unity property of shape functions at any point in the local domain
% Can be verified for xi, eta and zeta by varying (N_u(1 or 2 or 3,:))
actual_sum_N = sum(N_u(2,:));
expected_sum_N = 1;
try  
   cprintf('blue','Sum of shape functions test\n');
   verifyEqual(testCase, actual_sum_N, expected_sum_N, 'AbsTol', 1e-8);  
catch 
    disp('Sum of shape functions test failed');
end 
%% Test case: 2
%% Aim: To check vanishing property of shape function derivatives w.r.t  local coordinates Ni,xi & Ni,eta & Ni,zeta
% Can be verified for xi, eta and zeta by varying (N_u(1 or 2 or 3,:))
actual_sum_dN_u = sum(dN_u(2,:));
expected_sum_dN_u = 0;
try
    cprintf('blue','Sum of Ni,xi or Ni,eta or Ni,zeta test\n');
    verifyEqual(testCase, actual_sum_dN_u, expected_sum_dN_u, 'AbsTol', 1e-8);
catch
     disp('Sum of  shape function derivatives test');
end
%% Test case: 3
%% Aim: To check vanishing property of shape function derivatives w.r.t physical coordinates Ni,x & Ni,y & Ni,z
 actual_sum_dN_x = sum(B_f(2,:));
 expected_sum_dN_x = 0;
try
    cprintf('blue','Sum of Ni,x or Ni,y or Ni,z test\n');
    verifyEqual(testCase, actual_sum_dN_x, expected_sum_dN_x, 'AbsTol', 1e-8);
catch
     disp('Sum of  shape function derivatives test');
end
%% Test case: 4
%% Aim: To compute area/volume from determinant of jacobian
% Area/Volume
if problem == 1
   expected_area = 0.05;
   actual_area = area;
elseif problem == 3
   % Volume of cube of side 0.2, volume = 0.2^3
   expected_area = 0.008;
   actual_area = area;
else
   % Area of 2D element with size 0.2x0.2
   expected_area = 0.04;
   actual_area = area;
 end
try
     cprintf('blue','Length/Area/Volume of element test\n');
     verifyEqual(testCase, actual_area , expected_area, 'AbsTol', 1e-7);
catch
     disp('Area/volume of element test failed');
end
%% Test case: 5
%% Aim: To test the number of zero eigen values in element stiffness matrix
eigVals = real(eig(K_ele_mat + K_ele_geo));  relTol = 1e-12;
% use the maximum absolute real eigenvalue
refVal = max(abs(eigVals));
actual_zeros = sum(abs(eigVals) < relTol * refVal);
if problem == 1; expected_zeros = 3; elseif problem == 3; expected_zeros = 6; else; expected_zeros = 4; end
try
   cprintf('blue','Eigen value test\n');
   verifyEqual(testCase, actual_zeros , expected_zeros, 'AbsTol', 1e-8);
catch
    disp('Eigen value test failed');
end
%% Test case: 6
%% Aim: To calculate and test length/area while integrating external force on the boundary edge or surface
if problem == 2; len = 3 ; elseif problem == 3 ; len = 4; else; len = 2; end % No.of boundary nodes on element
if ~isempty(boundary_local_nodes) && length(boundary_local) == len 
if problem == 1
expected_ = 1;
actual_ = l_e;
elseif problem == 3
expected_ = 0.04;
actual_ = surface_area;
else
expected_ = 0.2;
actual_ = l_e;
end
try
   cprintf('blue','Length or area test on boundary\n');
   verifyEqual(testCase, actual_ , expected_, 'AbsTol', 1e-6);
catch    
    disp('Length or area test on boundary failed');
end
else
cprintf('blue','Length or area test on boundary\n');
fprintf('No boundary edge or surface\n');
fprintf('No test needed\n');
end
end
end

