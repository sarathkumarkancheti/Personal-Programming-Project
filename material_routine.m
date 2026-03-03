%% Material routine
%__________________________________________________________________________
% INPUT: Deformation gradient and its determinant, chemical potential of element, material parameters
%__________________________________________________________________________
% OUTPUT: Tangent modulus, cauchy stress and kirchoff stresses
%__________________________________________________________________________
function [D,C_stress,K_stress,S] = material_routine(F,J,cpe,G0,kB,NU,nu,X,ncoord) 


H1 = kB*NU/nu*(J*log(1-0.999/J)+1+X/J);            
H2 = kB*NU/nu*(J*log(1-0.999/J)+ J/(J-0.999)-X/J);
p_u = -cpe*kB*NU/nu +kB*NU/nu*(log(1-0.999/J)+1/J+X/J^2); % Volumetric term 

delta = eye(ncoord); % Kronecker delta
coeff1 = 2*G0-H1+cpe*kB*NU*J/nu;
coeff2 = H2-cpe*kB*NU*J/nu;
% Material tangent modulus
c = zeros(ncoord, ncoord, ncoord, ncoord);
for i = 1:ncoord
   for j = 1:ncoord
      for k = 1:ncoord
        for l = 1:ncoord
                    c(i,j,k,l) = coeff1*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))+coeff2*delta(i,j)*delta(k, l);
         end
      end
   end
end

C_stress = p_u*eye(ncoord)+ (G0/J)*(F*F.'-eye(ncoord));   % Cuachy stress
K_stress = J*C_stress;                                    % Kirchoff stress
S = J * inv(F)*C_stress * inv(F)';                        % Second Piola-Kirchhoff stress 
if ncoord == 2                                            % D matrix for 2D 
    D = [c(1,1,1,1), c(1,1,2,2),          0;
         c(2,2,1,1), c(2,2,2,2),          0;
                 0,           0,  c(1,2,1,2)];
elseif ncoord == 3                                        % D matrix for 3D
        D = [c(1,1,1,1), c(1,1,2,2), c(1,1,3,3), 0,          0,          0;
             c(2,2,1,1), c(2,2,2,2), c(2,2,3,3), 0,          0,          0;
             c(3,3,1,1), c(3,3,2,2), c(3,3,3,3), 0,          0,          0;
             0,          0,          0,          c(1,2,1,2), 0,          0;
             0,          0,          0,          0,          c(1,3,1,3), 0;
             0,          0,          0,          0,          0,          c(2,3,2,3)];
end
end