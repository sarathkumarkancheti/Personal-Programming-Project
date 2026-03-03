%% To calculate analytical solution for 1D problem(1)
%__________________________________________________________________________
% INPUT: Final time tend, number of time steps
%__________________________________________________________________________
%OUTPUT: Analytical stretch, chemical potential and stress
%__________________________________________________________________________
function [analytic_stretch,analytic_cp,analytic_sigma] = analytical_solution(tend,nt)
    H0 = 1;
    % Spatial and time mesh
    y = linspace(0, H0, 20);
    t = linspace(0, tend, nt);
    % Solve the PDE
    sol = pdepe(0, @pdefun, @icfun, @bcfun, y, t);
    % Extract the solution
    analytic_stretch = sol(:,:,1).';
    
    %Calculate chemical potential
    nu = 1.7e-28;  X = 0.2; G0 =1e7; kB = 1.3806488e-23;  theta = 298;
    
    analytic_cp = (G0*nu.*analytic_stretch - G0*nu*1./analytic_stretch+...
                        kB*theta*(log(1-0.999./analytic_stretch)+1./analytic_stretch + X./analytic_stretch.^2))/(kB*theta);
     
     
   analytic_sigma = -kB*theta/nu.*analytic_cp + kB*theta/nu*(log(1-0.999./analytic_stretch)+1./analytic_stretch + X./analytic_stretch.^2);
       
function [c, f, s] = pdefun(y, t, lambda, DlambdaDx)
    nu = 1.7e-28; 
     X = 0.2;
    G0 =1e7; 
    kB = 1.3806488e-23; 
    theta = 298;
    m_d = 0.5/(kB*nu*theta);
    % Define the PDE coefficients
    c = 1 / nu; % Coefficient
    f = m_d/lambda * (G0*nu*(1+1/lambda^2) + kB*theta*(1/(lambda*(lambda-0.999))-1/lambda^2-2*X/lambda^3))*DlambdaDx;
    s = 0;      % Source term
end
function lambda0 = icfun(y)
    % Initial condition
    if y < 1 %H0
        lambda0 = 1;
    else
        lambda0 = 1.498;
    end
end
function [pl, ql, pr, qr] = bcfun(yl, lambda_l, yr, lambda_r, t)
    % Boundary conditions
    pl = 0;
    ql = 1;
    pr = lambda_r-1.498 ;
    qr = 0;
end
end
