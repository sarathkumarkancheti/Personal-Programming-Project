 %% Function to calculate B matrices
%__________________________________________________________________________
% INPUT: Inverse jacobian, number of displacement and chemical potential
%        nodes for an element, number of coordinates, problem number
%__________________________________________________________________________
% OUTPUT: Gradient matrices in current configuration and G_grad in
%         reference configuration for defromation gradient
function [G,G_grad,B_u,B_cp,B_p,B_f] = B_matrices(invjacobian,invjacob_current,invjacob_current_c,dN_u,dN_cp,ncoord,nnode_ele_m,nnode_ele_c,problem) 
   %Multiplying with inverse jacobian for dN/dx, dN/dy
   %B_u, G, B_p and B_cp are evaluated in current configuration
   % G_grad for deformation gradient is evaluated in reference
   % configuration
   B = invjacob_current*dN_u;
   B_f = invjacobian*dN_u;
   B_c = invjacob_current_c*dN_cp;
   if problem == 1
      G = zeros(4,ncoord*nnode_ele_m);
      G_grad = zeros(4,ncoord*nnode_ele_m);
      B_u = zeros(3,ncoord*nnode_ele_m);
      B_cp = zeros(2,nnode_ele_c);
      B_p = zeros(1,ncoord*nnode_ele_m);
      for a = 1:nnode_ele_m
          %To calculate geometric stiffness K_geo
          G(1,2*a-1) = B(1,a); % Ni,x
          G(2,2*a-1) = B(2,a); % Ni,y
          G(3,2*a)   = B(1,a); % Ni,x
          G(4,2*a)   = B(2,a); % Ni,y
          %To calculate material stiffness K_mat
          B_u(1,2*a-1) = B(1,a);
          B_u(3,2*a-1) = B(2,a);
          B_u(2,2*a)   = B(2,a);
          B_u(3,2*a)   = B(1,a);
          %To calculate chemical stiffness H 
          B_cp(1,a) = B_c(1,a);
          B_cp(2,a) = B_c(2,a);
          %To calculate coupling stiffness P 
          B_p(2*a-1) = B(1,a);
          B_p(2*a) = B(2,a);
          G_grad(1,2*a-1) = B_f(1,a);
          G_grad(2,2*a-1) = B_f(2,a); 
          G_grad(3,2*a)   = B_f(1,a); 
          G_grad(4,2*a)   = B_f(2,a); 
       end
      
  elseif problem == 5 || problem == 2 
      G = zeros(4,ncoord*nnode_ele_m);
      G_grad = zeros(4,ncoord*nnode_ele_m);
      B_u = zeros(3,ncoord*nnode_ele_m);
      B_cp = zeros(2,nnode_ele_c);
      B_p = zeros(1,ncoord*nnode_ele_m);
      for a = 1:nnode_ele_m
          %To calculate geometric stiffness K_geo
          G(1,2*a-1) = B(1,a); % Ni,x
          G(2,2*a-1) = B(2,a); % Ni,y
          G(3,2*a)   = B(1,a); % Ni,x
          G(4,2*a)   = B(2,a); % Ni,y
          %To calculate material stiffness K_mat
          B_u(1,2*a-1) = B(1,a);
          B_u(3,2*a-1) = B(2,a);
          B_u(2,2*a)   = B(2,a);
          B_u(3,2*a)   = B(1,a);
          %To calculate chemical stiffness H 
          B_cp(1,a) = B_c(1,a);
          B_cp(2,a) = B_c(2,a);
          %To calculate coupling stiffness P 
          B_p(2*a-1) = B(1,a);
          B_p(2*a) = B(2,a);
          G_grad(1,2*a-1) = B_f(1,a); 
          G_grad(2,2*a-1) = B_f(2,a);
          G_grad(3,2*a)   = B_f(1,a); 
          G_grad(4,2*a)   = B_f(2,a); 
       end
   elseif problem == 3
      G = zeros(9,ncoord*nnode_ele_m);
      G_grad = zeros(9,24);
      B_u = zeros(6,ncoord*nnode_ele_m);
      B_cp = zeros(3,nnode_ele_m);
      B_p = zeros(1,ncoord*nnode_ele_m);
      for a = 1:nnode_ele_m
          %To calculate geometric stiffness K_geo
          G(1, 3*a-2) = B(1,a);  % u_x,x
          G(2, 3*a-2) = B(2,a);  % u_x,y
          G(3, 3*a-2) = B(3,a);  % u_x,z
          G(4, 3*a-1) = B(1,a);  % u_y,x
          G(5, 3*a-1) = B(2,a);  % u_y,y
          G(6, 3*a-1) = B(3,a);  % u_y,z
          G(7, 3*a)   = B(1,a);  % u_z,x
          G(8, 3*a)   = B(2,a);  % u_z,y
          G(9, 3*a)   = B(3,a);  % u_z,z
          %To calculate material stiffness K_mat
          B_u(1, 3*a-2) = B(1,a); 
          B_u(2, 3*a-1) = B(2,a);  
          B_u(3, 3*a)   = B(3,a);  
          B_u(4, 3*a-2) = B(2,a); 
          B_u(4, 3*a-1) = B(1,a);  
          B_u(5, 3*a-2) = B(3,a);  
          B_u(5, 3*a)   = B(1,a);  
          B_u(6, 3*a-1) = B(3,a);  
          B_u(6, 3*a)   = B(2,a); 
          %To calculate  chemical stiffness H 
          B_cp(1,a) = B_c(1,a); % Ni,x
          B_cp(2,a) = B_c(2,a); % Ni,y
          B_cp(3,a) = B_c(3,a); % Ni,z
          %To calculate coupling stiffness P 
          B_p(1,3*a-2) = B(1,a); % Ni,x
          B_p(1,3*a-1) = B(2,a); % Ni,y
          B_p(1,3*a)   = B(3,a); % Ni,z   
          G_grad(1, 3*a-2) = B_f(1,a);  % u_x,x
          G_grad(2, 3*a-2) = B_f(2,a);  % u_x,y
          G_grad(3, 3*a-2) = B_f(3,a);  % u_x,z
          G_grad(4, 3*a-1) = B_f(1,a);  % u_y,x
          G_grad(5, 3*a-1) = B_f(2,a);  % u_y,y
          G_grad(6, 3*a-1) = B_f(3,a);  % u_y,z
          G_grad(7, 3*a)   = B_f(1,a);  % u_z,x
          G_grad(8, 3*a)   = B_f(2,a);  % u_z,y
          G_grad(9, 3*a)   = B_f(3,a);  % u_z,z
      end
   elseif   problem == 4
      G = zeros(4,ncoord*nnode_ele_m);
      G_grad = zeros(4,ncoord*nnode_ele_m);
      B_u = zeros(3,ncoord*nnode_ele_m);
      B_cp = zeros(2,nnode_ele_c);
      B_p = zeros(1,ncoord*nnode_ele_m);
      for a = 1:nnode_ele_m
          %To calculate geometric stiffness K_geo
          G(1,2*a-1) = B(1,a); % Ni,x
          G(2,2*a-1) = B(2,a); % Ni,y
          G(3,2*a)   = B(1,a); % Ni,x
          G(4,2*a)   = B(2,a); % Ni,y
          %To calculate material stiffness K_mat
          B_u(1,2*a-1) = B(1,a);
          B_u(3,2*a-1) = B(2,a);
          B_u(2,2*a)   = B(2,a);
          B_u(3,2*a)   = B(1,a); 
          %To calculate coupling stiffness P 
          B_p(2*a-1) = B(1,a);
          B_p(2*a) = B(2,a);
          G_grad(1,2*a-1) = B_f(1,a);
          G_grad(2,2*a-1) = B_f(2,a); 
          G_grad(3,2*a)   = B_f(1,a); 
          G_grad(4,2*a)   = B_f(2,a); 
      end
       for b = 1:nnode_ele_c
           %To calculate chemical stiffness H 
          B_cp(1,b) = B_c(1,b);
          B_cp(2,b) = B_c(2,b);
       end
   end
end
