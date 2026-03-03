%% Material parameters
%__________________________________________________________________________
% INPUT: Nominal diffusion constant, number of coordinates(2 for 2D, 3 for 3D)
%__________________________________________________________________________
% OUTPUT: G0,kB,NU,nu,X,md,mc,m
%__________________________________________________________________________
function [G0,kB,NU,nu,X,md,mc,m] = material_parameters(m_d,ncoord)
  G0 = 1e7;              % Shear modulus [Pa]
  kB = 1.3806488e-23;    % Boltzmann constant [J/K]
  NU = 298;              % Absolute temperature [K]
  nu = 1.7e-28;          % Volume of solvent molecule [m^3]
  X = 0.2;               % Interaction dimensionless parameter
  m_c = 0.1;             % Material convection constant
  md = m_d/(kB*NU*nu);   % Diffusion coefficient
  mc = m_c/(kB*NU*nu);   % Convection coefficient
  m = md*eye(ncoord);    % Mobility tensor
end