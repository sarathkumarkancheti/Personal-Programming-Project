                                    %% Shape function derivatives
%__________________________________________________________________________
%INPUT: Local coordinates(xi,eta,zeta,problem)
%__________________________________________________________________________
%OUTPUT:[dN_u] for displacements and [dN_cp] for chemical potential
%__________________________________________________________________________
function [dN_u,dN_cp] = shape_function_derivatives(xi,eta,zeta,problem)
     if problem ==1
     % 4 nodes for both displacements and chemical potential
     % Derivatives w.r.t xi   
     dN1_xi = -(1/4)*(1-eta);
     dN2_xi = (1/4)*(1-eta);
     dN3_xi = (1/4)*(1+eta);
     dN4_xi = -(1/4)*(1+eta); 
     % Derivatives w.r.t eta
     dN1_eta = -(1/4)*(1-xi);
     dN2_eta = -(1/4)*(1+xi);
     dN3_eta = (1/4)*(1+xi);
     dN4_eta = (1/4)*(1-xi);
     
     dN_xi = [dN1_xi,dN2_xi,dN3_xi,dN4_xi];       %[N],ξ
     dN_eta = [dN1_eta,dN2_eta,dN3_eta,dN4_eta];  %[N],η
     dN_u = [dN_xi;dN_eta];
     dN_cp = [dN_xi;dN_eta];
     
     elseif problem == 2 || problem == 5
     % 8 nodes for both displacements and chemical potential
     % Derivatives w.r.t xi
     dN1_xi = -0.25*(-1+eta)*(2*xi+eta);
     dN2_xi =  0.25*(-1+eta)*(eta-2*xi);
     dN3_xi =  0.25*(eta+1)*(2*xi+eta);
     dN4_xi =  0.25*(eta+1)*(2*xi-eta);
     dN5_xi =  xi*(-1+eta);
     dN6_xi =  0.5*(1-eta^2);
     dN7_xi =  -xi*(1+eta); 
     dN8_xi =  -0.5*(1-eta^2);
     % Derivatives w.r.t eta
     dN1_eta = -0.25*(xi-1)*(xi+2*eta);
     dN2_eta =  0.25*(1+xi)*(2*eta-xi);
     dN3_eta =  0.25*(xi+1)*(xi+2*eta);
     dN4_eta = -0.25*(-1+xi)*(2*eta-xi);
     dN5_eta =  0.5*(xi^2-1);
     dN6_eta =  -eta*(xi+1);
     dN7_eta =  0.5*(1-xi^2);
     dN8_eta =  -eta*(1-xi);
   
     dN_xi = [dN1_xi,dN2_xi,dN3_xi,dN4_xi,dN5_xi,dN6_xi,dN7_xi,dN8_xi];          %[N],ξ
     dN_eta = [dN1_eta,dN2_eta,dN3_eta,dN4_eta,dN5_eta,dN6_eta,dN7_eta,dN8_eta]; %[N],η
     dN_u = [dN_xi;dN_eta];
     dN_cp = [dN_xi;dN_eta];
    
     elseif problem ==3
     % 8 nodes for both displacements and chemical potential
     % Derivatives w.r.t xi
     dN1_xi = -1/8*(1-eta)*(1-zeta);
     dN2_xi =  1/8*(1-eta)*(1-zeta);
     dN3_xi =  1/8*(1+eta)*(1-zeta);
     dN4_xi = -1/8*(1+eta)*(1-zeta);
     dN5_xi = -1/8*(1-eta)*(1+zeta);
     dN6_xi =  1/8*(1-eta)*(1+zeta);
     dN7_xi =  1/8*(1+eta)*(1+zeta);
     dN8_xi = -1/8*(1+eta)*(1+zeta);
     
     % Derivatives w.r.t eta
     dN1_eta= -1/8*(1-xi)*(1-zeta);
     dN2_eta= -1/8*(1+xi)*(1-zeta);
     dN3_eta= 1/8*(1+xi)*(1-zeta);
     dN4_eta= 1/8*(1-xi)*(1-zeta);
     dN5_eta= -1/8*(1-xi)*(1+zeta);
     dN6_eta= -1/8*(1+xi)*(1+zeta);
     dN7_eta= 1/8*(1+xi)*(1+zeta);
     dN8_eta= 1/8*(1-xi)*(1+zeta);
     
     % Derivatives with respect to zeta
     dN1_zeta= -1/8*(1-xi)*(1-eta);
     dN2_zeta= -1/8*(1+xi)*(1-eta);
     dN3_zeta= -1/8*(1+xi)*(1+eta);
     dN4_zeta= -1/8*(1-xi)*(1+eta);
     dN5_zeta= 1/8*(1-xi)*(1-eta);
     dN6_zeta= 1/8*(1+xi)*(1-eta);
     dN7_zeta= 1/8*(1+xi)*(1+eta);
     dN8_zeta= 1/8*(1-xi)*(1+eta);
     
     dN_xi = [dN1_xi,dN2_xi,dN3_xi,dN4_xi,dN5_xi,dN6_xi,dN7_xi,dN8_xi];                   %[N],ξ
     dN_eta = [dN1_eta,dN2_eta,dN3_eta,dN4_eta,dN5_eta,dN6_eta,dN7_eta,dN8_eta];          %[N],η
     dN_zeta = [dN1_zeta,dN2_zeta,dN3_zeta,dN4_zeta,dN5_zeta,dN6_zeta,dN7_zeta,dN8_zeta]; %[N],ζ
     
     dN_u = [dN_xi;dN_eta;dN_zeta];
     dN_cp = [dN_xi;dN_eta;dN_zeta];
     
     elseif problem == 4 
     % 8 nodes for displacements and 4 nodes for chemical potential
     % Derivatives w.r.t xi   
     dN1_xi = -(1/4)*(1-eta);
     dN2_xi = (1/4)*(1-eta);
     dN3_xi = (1/4)*(1+eta);
     dN4_xi = -(1/4)*(1+eta); 
     % Derivatives w.r.t eta
     dN1_eta = -(1/4)*(1-xi);
     dN2_eta = -(1/4)*(1+xi);
     dN3_eta = (1/4)*(1+xi);
     dN4_eta = (1/4)*(1-xi);
     dN_xi = [dN1_xi,dN2_xi,dN3_xi,dN4_xi];      %[N],ξ
     dN_eta = [dN1_eta,dN2_eta,dN3_eta,dN4_eta]; %[N],η
     dN_cp = [dN_xi;dN_eta];
     
     % Derivatives w.r.t xi
     dN1m_xi = -0.25*(-1+eta)*(2*xi+eta);
     dN2m_xi =  0.25*(-1+eta)*(eta-2*xi);
     dN3m_xi =  0.25*(eta+1)*(2*xi+eta);
     dN4m_xi =  0.25*(eta+1)*(2*xi-eta);
     dN5m_xi =  xi*(-1+eta);
     dN6m_xi =  0.5*(1-eta^2);
     dN7m_xi =  -xi*(1+eta); 
     dN8m_xi =  -0.5*(1-eta^2);
     
     % Derivatives w.r.t eta
     dN1m_eta = -0.25*(xi-1)*(xi+2*eta);
     dN2m_eta =  0.25*(1+xi)*(2*eta-xi);
     dN3m_eta =  0.25*(xi+1)*(xi+2*eta);
     dN4m_eta = -0.25*(-1+xi)*(2*eta-xi);
     dN5m_eta =  0.5*(xi^2-1);
     dN6m_eta =  -eta*(xi+1);
     dN7m_eta =  0.5*(1-xi^2);
     dN8m_eta =  -eta*(1-xi);
   
     dNm_xi = [dN1m_xi,dN2m_xi,dN3m_xi,dN4m_xi,dN5m_xi,dN6m_xi,dN7m_xi,dN8m_xi];          %[N],ξ
     dNm_eta = [dN1m_eta,dN2m_eta,dN3m_eta,dN4m_eta,dN5m_eta,dN6m_eta,dN7m_eta,dN8m_eta]; %[N],η
     dN_u = [dNm_xi;dNm_eta];
     end
     
end
