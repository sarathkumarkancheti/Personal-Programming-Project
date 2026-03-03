%% Shape functions
%__________________________________________________________________________
%INPUT: Local coordinates(xi,eta,zeta,problem)
%__________________________________________________________________________
%OUTPUT:[N_u] for displacements and [N_cp] for chemical potential
%__________________________________________________________________________
function [N_u,N_cp] = shape_functions(xi,eta,zeta,problem)
     if problem == 1
     % For Q4μ4 element
     N1 = (1/4)*(1-xi)*(1-eta);
     N2 = (1/4)*(1+xi)*(1-eta);
     N3 = (1/4)*(1+xi)*(1+eta);
     N4 = (1/4)*(1-xi)*(1+eta);
     % Shape functions for displacements
     N_u = [[N1,0,N2,0,N3,0,N4,0];...
            [0,N1,0,N2,0,N3,0,N4]];
     % Shape functions for chemical potential
     N_cp = [N1,N2,N3,N4];
        
     elseif problem ==2 || problem == 5
     % For Q8μ8 element
     N1 = -0.25*(1-xi)*(1-eta)*(1+xi+eta);
     N2 = -0.25*(1+xi)*(1-eta)*(1-xi+eta);
     N3 = -0.25*(1+xi)*(1+eta)*(1-xi-eta);
     N4 = -0.25*(1-xi)*(1+eta)*(1+xi-eta);
     N5 = 0.5*(1-xi^2)*(1-eta);
     N6 = 0.5*(1+xi)*(1-eta^2);
     N7 = 0.5*(1-xi^2)*(1+eta);
     N8 = 0.5*(1-xi)*(1-eta^2);   
     % Shape functions for displacements
     N_u = [[N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8,0];
            [0,N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8]];   
     % Shape functions for chemical potential
     N_cp = [N1,N2,N3,N4,N5,N6,N7,N8];
   
     elseif problem ==3
     % For H8μ8 element
     N1 = 1/8*(1-xi)*(1-eta)*(1-zeta);
     N2 = 1/8*(1+xi)*(1-eta)*(1-zeta);
     N3 = 1/8*(1+xi)*(1+eta)*(1-zeta);
     N4 = 1/8*(1-xi)*(1+eta)*(1-zeta);
     N5 = 1/8*(1-xi)*(1-eta)*(1+zeta);
     N6 = 1/8*(1+xi)*(1-eta)*(1+zeta);
     N7 = 1/8*(1+xi)*(1+eta)*(1+zeta);
     N8 = 1/8*(1-xi)*(1+eta)*(1+zeta);
     % Shape functions for displacements
     N_u = [[N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0,0];...
           [0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0];...
           [0,0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8]];
     % Shape functions for chemical potential
     N_cp = [N1,N2,N3,N4,N5,N6,N7,N8];  
     
     elseif problem == 4  
     % For Q8μ4 element
     N1 = -0.25*(1-xi)*(1-eta)*(1+xi+eta);
     N2 = -0.25*(1+xi)*(1-eta)*(1-xi+eta);
     N3 = -0.25*(1+xi)*(1+eta)*(1-xi-eta);
     N4 = -0.25*(1-xi)*(1+eta)*(1+xi-eta);
     N5 = 0.5*(1-xi^2)*(1-eta);
     N6 = 0.5*(1+xi)*(1-eta^2);
     N7 = 0.5*(1-xi^2)*(1+eta);
     N8 = 0.5*(1-xi)*(1-eta^2); 
     % Shape functions for 8 displacement nodes
     N_u = [[N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8,0];
            [0,N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8]];
   
     N1 = (1/4)*(1-xi)*(1-eta);
     N2 = (1/4)*(1+xi)*(1-eta);
     N3 = (1/4)*(1+xi)*(1+eta);
     N4 = (1/4)*(1-xi)*(1+eta);
    % Shape functions for 4 chemical potential nodes
     N_cp = [N1,N2,N3,N4];
     end
     
end

    