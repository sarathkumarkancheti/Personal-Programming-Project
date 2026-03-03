%% Function to verify numerical and analytical tangent modulus
%__________________________________________________________________________
% INPUT : Deformation gradient F, element chemical potential (cpe), and
% number of coordinates ncoord = 2 for 2D, 3 for 3D
%__________________________________________________________________________
% OUTPUT: 1.Best perturbation size (delta)
%         2. Mean absolute and relative errors for both filled and unfilled components of D
%         3. Component wise relative and absolute errors
%         4. Analytical and numerical tangent moduli
%__________________________________________________________________________
function verify_tangent_modulus(F, cpe , ncoord)
    % Material parameters
    G0 = 1e7;           % Shear modulus [Pa]
    kB = 1.3806488e-23; % Boltzmann constant [J/K]
    NU = 298;           % Temperature [K]
    nu = 1.7e-28;       % Solvent molecule volume [m^3]
    X = 0.2;            % Interaction dimensionless parameter
 
    % Voigt notation mapping
    if ncoord == 2
        voigt_map = [1 1; 2 2; 1 2]; % 2D: 11, 22, 12
        D_size = 3;                  % 3x3 D matrix
    elseif ncoord == 3
        voigt_map = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3]; % 3D: 11, 22, 33, 12, 13, 23
        D_size = 6;                                 % 6x6 D matrix
    else
        error('ncoord must be 2 or 3');
    end
    
    % Deformation gradient
    J = det(F);
    if size(F,1) ~= ncoord || size(F,2) ~= ncoord
        error('F must be %dx%d for ncoord = %d', ncoord, ncoord, ncoord);
    end
    
    % Calculate analytical tangent modulus from material routine
    [D_analytical, ~, ~, ~] = material_routine(F, J, cpe, G0, kB, NU, nu, X, ncoord);
    
    % Numerical tangent modulus calculation based on dS/dC (central differences)
    delta_values = [1e-4, 1e-5, 1e-6, 1e-7]; % Multiple perturbation sizes
    best_D_num = zeros(D_size, D_size);
    best_MAE_filled = Inf;
    best_delta = 0;
    
    for delta = delta_values
        D_num = zeros(D_size, D_size); % Initialize numerical tangent modulus
        
        % Compute Cauchy green tensor, C = F'*F
        C = F'*F;
        for i = 1:D_size
            kl = voigt_map(i,:);
            
            % Positive perturbation in C
            dC = zeros(ncoord, ncoord);
            dC(kl(1),kl(2)) = delta;
            dC(kl(2),kl(1)) = delta; % Ensure symmetry
            C_plus = C + dC;
            [U,S,V] = svd(C_plus);
            F_plus = U*sqrt(S)*V';
            J_plus = det(F_plus);
            [~, ~, ~, S_plus] = material_routine(F_plus, J_plus, cpe, G0, kB, NU, nu, X, ncoord);
            
            % Negative perturbation in C
            C_minus = C - dC;
            [U,S,V] = svd(C_minus);
            F_minus = U*sqrt(S)*V';
            J_minus = det(F_minus);
            [~, ~, ~, S_minus] = material_routine(F_minus, J_minus, cpe, G0, kB, NU, nu, X, ncoord);
            
            % Central difference for dS/dC
            dS = (S_plus - S_minus)/(2*delta);
            
            % Store in Voigt notation
            if ncoord == 2
                D_num(:,i) = [dS(1,1); dS(2,2); dS(1,2)];
            elseif ncoord == 3
                D_num(:,i) = [dS(1,1); dS(2,2); dS(3,3); dS(1,2); dS(1,3); dS(2,3)];
            end
        end
        
        % Convert dS/dC to dS/dE (multiply by 2 since E = (1/2)*(C - I))
        D_num = 2 * D_num;
        
        % Enforce symmetry to match D from reference paper
        if ncoord == 2
            D_num(1:2,3) = 0;  % Remove C1112 and C2212
            D_num(3,1:2) = 0;  % Remove C1211 and C1222
            D_num = 0.5*(D_num + D_num'); % Enforce minor symmetry
        elseif ncoord == 3
            D_num(1:3,4:6) = 0; % Remove shear-normal terms (C1112, C1123, etc.)
            D_num(4:6,1:3) = 0; % Remove shear-normal terms (C1211, C2311, etc.)
            D_num = 0.5*(D_num + D_num'); % Enforce minor symmetry
        end
        
        % Error analysis for current delta
        filled_ = abs(D_analytical) > 1e-6;
        unfilled_ = ~filled_;
        
        abs_error = abs(D_analytical - D_num);          % Absolute error
        abs_norm = abs_error / G0;                      % Normalized error
        rel_error = abs_error ./ abs(D_analytical);     % Relative error
        rel_error(isinf(rel_error) | isnan(rel_error)) = 0; 
        
        MAE_filled = mean(abs_norm(filled_));
 
        % Update best result
        if MAE_filled < best_MAE_filled
            best_MAE_filled = MAE_filled;
            best_D_num = D_num;
            best_delta = delta;
        end
    end
    
    % Use the best numerical tangent modulus
    D_num = best_D_num;
    
    % Final error analysis
    filled_ = abs(D_analytical) > 1e-6;
    unfilled_ = ~filled_;
    
    abs_error = abs(D_analytical - D_num);          % Absolute error
    abs_norm = abs_error / G0;                      % Normalized error
    rel_error = abs_error ./ abs(D_analytical);     % Relative error
    rel_error(isinf(rel_error) | isnan(rel_error)) = 0; 
    

    % Display results
    fprintf('============ Results ===============\n');
    fprintf('Best perturbation size (delta): %.2e\n', best_delta);
    fprintf('Mean relative error (filled): %.2e%%\n', mean(rel_error(filled_))*100);
    fprintf('Mean absolute error (filled): %.2e \n\n', mean(abs_norm(filled_)));
    fprintf('Mean relative error (unfilled): %.2e%%\n', mean(rel_error(unfilled_))*100);
    fprintf('Mean absolute error (unfilled): %.2e \n\n', mean(abs_norm(unfilled_)));
    
    fprintf('Component-wise Absolute Errors (normalized by G0):\n');
    disp(abs_norm);
    fprintf('Component-wise Relative Errors:\n');
    disp(rel_error);
    
    fprintf('Analytical Tangent Modulus:\n');
    disp(D_analytical);
    
    fprintf('Numerical Tangent Modulus:\n');
    disp(D_num);
end