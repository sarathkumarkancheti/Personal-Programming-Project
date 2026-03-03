%% Jacobian matrix [dx/dξ, dx/dη ; dy/dξ, dy/dη]
%__________________________________________________________________________
%INPUT:Shape function derivatives [dN,ξ], [dN,η], and element coordinates 
%__________________________________________________________________________
%OUTPUT:Jacobian matrix, inverse jacobian matrix and determinant of
%       jacobian for reference and current configurations
%__________________________________________________________________________
function [Jacobian,invjacobian,detjacobian,invjacob_current,detjacobian_c,invjacob_current_c] = jacobian(nnode_ele_c,dN_u,dN_cp,ele_coord,element_displacement)
     % Current cooedinates of element
     el_coord_current = ele_coord +  reshape(element_displacement,size(ele_coord,2),size(ele_coord,1)).';
     Jacobian = dN_u*(ele_coord);                  % Reference configuration jacobian
     detjacobian = det(Jacobian);
     invjacobian = inv(Jacobian);
     Jacobian_current = dN_u*el_coord_current;     % Current configuration jacobian
     invjacob_current = inv(Jacobian_current);
     
     Jacobian_c = dN_cp*(ele_coord(1:nnode_ele_c,:));
     detjacobian_c = det(Jacobian_c);
     %invjacobian_c = inv(Jacobian_c);
     Jacobian_current_c = dN_cp*el_coord_current(1:nnode_ele_c,:);
     invjacob_current_c = inv(Jacobian_current_c);
     size(Jacobian);
end
