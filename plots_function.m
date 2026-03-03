function plots_function(problem)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots chemical potential, stretch and stress for both numerical and
% analytical solution
if problem == 1    
% Read analytical solution data for tend = 5, number of time steps = 1001
% with time increment Δt = 0.005
%[analytic_stretch,analytic_cp,analytic_sigma] = analytical_solution(5,1001);
% Load chemical potential
mu = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\chemical_potential.dat",'\t',1,0);
% Load stretch  
stretch = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\stretch.dat",'\t',1,0); 
% Load stress 
sigma_xx = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\stress.dat",'\t',1,0);  

mu_ref = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\mu_data.csv"); % Extracted using WebPlotDigitizer
lambda_ref = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\lambda_data.csv"); % Extracted using WebPlotDigitizer
stress_ref = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\stress_data.csv"); % Extracted using WebPlotDigitizer
x_ref1 = mu_ref(:,1); % mu_ref spatial points from WebPlotDigitizer
y_ref1 = mu_ref(:,2); % mu_ref at different times from WebPlotDigitizer
x_lam1 = lambda_ref(:,1); % stretch_ref spatial points from WebPlotDigitizer
y_lam1 = lambda_ref(:,2); % stretch_ref at different times from WebPlotDigitizer
x_stress = stress_ref(:,1); % stress_ref spatial points from WebPlotDigitizer
y_stress = stress_ref(:,2); % stress_ref at different times from WebPlotDigitizer

%======================== Plot for chemical potential =====================
%__________________________________________________________________________
all_time_steps = linspace(0, 5, 1001); 
selected_times = [0,0.01,0.04,0.19,0.39,1.0,5];  % Times to plot
[~, indices] = min(abs(all_time_steps'-selected_times),[],1); % Find nearest indices
X2 = linspace(0, 1, size(mu, 1));

figure;
hold on;
colors = lines(length(selected_times)); 
h = gobjects(length(selected_times), 1);  % Preallocate line handles

for i = 1:length(selected_times)
    idx = indices(i);
    % Numerical Solution
   smoothed_mu = smooth(mu(:, idx), 1);
    h(i) = plot(X2, smoothed_mu, 'Color', colors(i, :), 'LineWidth', 1.5);
    Y(:,i) = smoothed_mu;
end
ha = plot(x_ref1, y_ref1, 'ko', 'MarkerSize', 5, 'DisplayName', 'Reference');
xlabel('$Y/H_0$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\hat{\mu}$', 'Interpreter', 'latex', 'FontSize', 10);
legend([h; ha], [arrayfun(@(t) sprintf('t = %.2f', t), selected_times, 'UniformOutput', false), {'Reference'}], ...
       'Location', 'southwest', 'FontSize', 8);
annotation('arrow', [0.9 0.3], [0.2 0.9],'LineWidth', 1.5);
set(gca,'FontSize',10);
ylim([-6, 1]);
text(0.7, 0.2, 'Steady state', 'FontSize', 10, 'FontWeight', 'bold');
sgtitle('Nominal chemical potential','Fontsize',10);
grid on; box on;
hold off;

%========================== Plot for stretch ==============================
%__________________________________________________________________________
X2 = linspace(0, 1, size(stretch, 1));
figure;
hold on;
gap_center = 0.5;        % X-location of the gap
gap_width = 0.05;         % Width of the gap
% Plot reference data
ha = plot(x_lam1, y_lam1, 'ko', 'MarkerSize', 5, 'DisplayName', 'Reference');
for i = 1:length(selected_times)
idx = indices(i);
xdata = X2;
ydata = stretch(:, idx);
left = xdata < (gap_center - gap_width/2);
right = xdata > (gap_center + gap_width/2);
 % Plot left and right parts of the line separately
h(i) = plot(xdata(left), ydata(left), 'Color', colors(i,:), 'LineWidth', 1.5);
            plot(xdata(right), ydata(right), 'Color', colors(i,:), 'LineWidth', 1.5);
Y_lam(:,i) = stretch(:,idx);
% Add text label at the center of the gap
y_gap = interp1(xdata, ydata, gap_center); % estimate y at gap center
box_width = 0.15;
box_height = 0.03;
% Draw white rectangle behind the text (in data units)
rectangle('Position', [gap_center - box_width/2, y_gap - box_height/2, box_width, box_height], ...
          'FaceColor', 'white', 'EdgeColor', 'none');
% Draw the text on top of the rectangle
text(gap_center, y_gap, sprintf('t = %.2f', selected_times(i)), ...
     'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
xlabel('$Y/H_0$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 10);
legend([h(1), ha], {'Numerical', 'Reference'},'Location', 'southwest', 'FontSize', 8);
annotation('arrow', [0.9 0.3], [0.2 0.9],'LineWidth', 1.5);
text(0.7, 1.52, 'Steady state', 'FontSize', 10, 'FontWeight', 'bold');
sgtitle('Stretch in Y direction','Fontsize',10);
ylim([0.9, 1.6]);
grid on; box on; hold off;

%============================= Plot for stress ============================
%__________________________________________________________________________
X2 = linspace(0, 1, size(sigma_xx, 1));
figure;
hold on;

gap_center = 0.5;
gap_width = 0.15;
gap_left = gap_center - gap_width/2;
gap_right = gap_center + gap_width/2;

% Split reference markers
left_mask = x_stress < gap_left;
right_mask = x_stress > gap_right;

% Plot reference markers only outside the gap
ha1 = plot(x_stress(left_mask), y_stress(left_mask), 'ko', 'MarkerSize', 5, 'DisplayName', 'Reference');
ha2 = plot(x_stress(right_mask), y_stress(right_mask), 'ko', 'MarkerSize', 5, 'HandleVisibility', 'off');

% Plot lines with gaps and insert text
for i = 1:length(selected_times)
    idx = indices(i);  
    xdata = X2;
    ydata = sigma_xx(:, idx);
    Y_stress(:,i) = sigma_xx(:,idx);
    left = xdata < gap_left;
    right = xdata > gap_right;
    % Plot left and right segments
    h(i) = plot(xdata(left), ydata(left), 'Color', colors(i,:), 'LineWidth', 1.5);
            plot(xdata(right), ydata(right), 'Color', colors(i,:), 'LineWidth', 1.5);
    % Interpolated value at gap center
    y_gap = interp1(xdata, ydata, gap_center);
    box_width = 0.3;
    box_height = 0.5;
    rectangle('Position', [gap_center - box_width/2, y_gap - box_height/2, box_width, box_height], ...
              'FaceColor', 'white', 'EdgeColor', 'none');
    text(gap_center, y_gap, sprintf('t = %.2f', selected_times(i)), ...
         'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
xlabel('$Y/H_0$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\sigma_X$', 'Interpreter', 'latex', 'FontSize', 10);
ylim([-10e6, 1e6]);
legend([h(1), ha1], {'Numerical', 'Reference'}, 'Location', 'southwest', 'FontSize', 9);
annotation('arrow', [0.2 0.8], [0.9 0.18],'LineWidth', 1.5);
text(0.6, -8.5e6, 'Steady state', 'FontSize', 10, 'FontWeight', 'bold', 'Rotation', 0);
sgtitle('Evolution of stress','Fontsize',10);
grid on; box on;
hold off;

%========================== Error calculation==============================
times = [0.01,0.04,0.19,0.39,1.0];
folder = 'C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_1\ref_data';
for i = 1:5
    filename = sprintf('mu_%d.csv', i);      % Read chemical potential from WebPlotDigitizer
    filename1 = sprintf('lambda_%d.csv', i); % Read stretch from WebPlotDigitizer
    filename2 = sprintf('sigma_%d.csv', i);  % Read stress from WebPlotDigitizer
    filepath = fullfile(folder, filename);   
    filepath1 = fullfile(folder, filename1); 
    filepath2 = fullfile(folder, filename2);
    if isfile(filepath)
    data = csvread(filepath);
    data1 = csvread(filepath1);
    data2 = csvread(filepath2);
    x_ref = data(:,1);
    y_ref = data(:,2);
    cp_num = Y(:,i+1);
    y_mu_interp = interp1(X2,cp_num,x_ref,'linear');
    x_lam = data1(:,1);
    y_lam = data1(:,2);
    y_num = Y_lam(:,i+1);
    y_lam_interp = interp1(X2,y_num,x_lam,'linear');
    x_sig = data2(:,1);
    y_sig = data2(:,2);
    sig_num = Y_stress(:,i+1);
    y_sig_interp = interp1(X2,sig_num,x_sig,'linear');
    error_cp = y_ref - y_mu_interp;
    error_lam = y_lam - y_lam_interp;
    error_sig = (y_sig - y_sig_interp)/1e6;
    RMSE1 = sqrt(mean(error_cp.^2));
    MAE1 = mean(abs(error_cp));
    RMSE2 = sqrt(mean(error_lam.^2));
    MAE2 = mean(abs(error_lam));
    RMSE3 = sqrt(mean(error_sig.^2));
    MAE3 = mean(abs(error_sig));
   fprintf('Time: %f\n', times(i));
   fprintf('\x03BC --> RMSE: %.6f, MAE: %.6f\n', RMSE1, MAE1);
   fprintf('\x03BB --> RMSE: %.6f, MAE: %.6f\n', RMSE2, MAE2);
   fprintf('\x03C3 --> RMSE: %.6f, MAE: %.6f\n', RMSE3, MAE3);
    else
        fprintf('File not found: %s\n', filepath);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots y-displacements of points C,D and stress_xx in deformed
% configuration
elseif problem == 2
    
% Load the displacement vector
displacement = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\displacement_vector.dat", '\t', 1, 0);
% Load Gauss point stress 
stress_xx_gp = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\stress_xx.dat", '\t', 1, 0);
stress_all = reshape(stress_xx_gp, [25, 4, 1001]); %reshape into orginal size
error = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\error.dat",'\t',1,0);
% Load mesh data
[~, ~, ~, elementNodes, ~, ~, ~, ~, nodeCoords, ~] = inp_coordinates(2,2);
node_coords = [nodeCoords(:,1), nodeCoords(:,2)];
disp_ref2 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\y_c_data.csv"); % Extracted using WebPlotDigitizer
disp_ref1 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\y_d_data.csv"); % Extracted using WebPlotDigitizer

%=========================== Convergence plot =============================
%__________________________________________________________________________
t = 11;  % First 10 time steps
[~, total_time_steps] = size(error);
if total_time_steps < t
    error("Your 'error' matrix has fewer than 10 time steps.");
end
figure;
hold on;
for j = 1:t
    err_vals = error(:, j+1);
    % Remove zeros (assumes iterations stop after convergence)
    valid_idx = err_vals ~= 0;
    err_vals = err_vals(valid_idx);
    iter_count = numel(err_vals);
    h = plot(j * ones(1, iter_count), err_vals, 'k-o','MarkerFaceColor', 'k', 'MarkerSize', 7);
end
legend(h, 'Iterations', 'Location', 'northeast');
set(gca, 'YScale', 'log');
xlabel('Time steps');
ylabel('e_{disp}');
xlim([0 t+1]);
sgtitle('Convergence iterations','Fontsize',10);
grid on; box on;
%title('Error indicator evolution over first 10 time steps');


%================== Plot for y coordinates of points C and D ==============
%__________________________________________________________________________
x = linspace(0,5,1001);  %tend = 5, nt = 1001 with time increment Δt = 0.005
y1 = nodeCoords(31,2)+displacement(2*31,:); %y coordinate of D (node 31 in mesh)
y2 = nodeCoords(36,2)+displacement(2*36,:); %y coordinate of C (node 36 in mesh)
% Reference data
x_ref1 = disp_ref1(:,1); % D time points extracted from WebPlotDigitizer
y_ref1 = disp_ref1(:,2); % D displacement points extracted from WebPlotDigitizer
x_ref2 = disp_ref2(:,1); % C time points extracted from WebPlotDigitizer
y_ref2 = disp_ref2(:,2); % C displacement points extracted from WebPlotDigitizer
% Plot the y-coordinates
figure;
hold on;
% Numerical solution
plot(x, y1, '-r', 'LineWidth', 1.5, 'DisplayName', 'y_D');
plot(x, y2, '-b', 'LineWidth', 1.5, 'DisplayName', 'y_C');

% Reference solution(Markers)
[x_ref1_sorted, sort_idx1] = sort(x_ref1);
y_ref1_sorted = y_ref1(sort_idx1);
[x_ref2_sorted, sort_idx2] = sort(x_ref2);
y_ref2_sorted = y_ref2(sort_idx2);
plot(x_ref1_sorted, y_ref1_sorted, 'ko', 'MarkerSize', 5, 'DisplayName', 'Reference data');
plot(x_ref2_sorted, y_ref2_sorted, 'ko','MarkerSize', 5, 'HandleVisibility', 'off');

% Vertical lines
t1 = 1.00;
plot([t1 t1], [0.9 y1(round(t1/0.005)+1)], '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t1+0.05, 0.92, ['t = ' num2str(t1)], 'HorizontalAlignment', 'center', 'FontSize', 10);
plot(t1, y1(round(t1/0.005)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t1, y2(round(t1/0.005)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');

t2 = 1.94;
plot([t2 t2], [0.9 y1(round(t2/0.005)+1)], '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t2+0.05, 1.1, ['t = ' num2str(t2)], 'HorizontalAlignment', 'center', 'FontSize', 10);
plot(t2, y1(round(t2/0.005)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t2, y2(round(t2/0.005)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');

t3 = 5;
plot(t3, (y1(round(t3/0.005)+1)+y2(round(t3/0.005)+1))/2, 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');

t4 = 0;
plot(t4, (y1(round(t4/0.005)+1)+y2(round(t4/0.005)+1))/2, 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
hold off;
box on;
xlabel('Time','FontSize', 10);
ylabel('y Coordinates','FontSize', 10);
legend('Location', 'northeast', 'FontSize', 10);
sgtitle('Evolutions of y coordinates of points C and D','Fontsize',10);
xlim([0 5]); ylim([0.9 1.5]);

%========================= Error calculation ==============================
%__________________________________________________________________________
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_2\error_results_2.txt', 'w');
% Interpolate yC and yD for the reference time points 
y_1_interp = interp1(x,y1,x_ref1,'linear');
y_2_interp = interp1(x,y2,x_ref2,'linear');
% Calculate error
error1 = y_1_interp - y_ref1;
error2 = y_2_interp - y_ref2;
%Calculate Root mean squared error , Mean absolute error for yD and yC
RMSE1 = sqrt(mean(error1.^2));
MAE1 = mean(abs(error1));
RMSE2 = sqrt(mean(error2.^2));
MAE2 = mean(abs(error2));
RMSE1_percent = (RMSE1/mean(y_ref1))*100;
MAE1_percent  = (MAE1 /mean(y_ref1))*100;
RMSE2_percent = (RMSE2/ mean(y_ref2))*100;
MAE2_percent  = (MAE2 / mean(y_ref2))*100;
fprintf('y_C --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf('y_D --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
% Write to file
fprintf(fileID,'y_C --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf(fileID,'y_D --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
fclose(fileID);

%=========================== Plot for stress ==============================
%__________________________________________________________________________
% Gauss point locations for 2x2 integration
A = [ -1 -1;1 -1;1  1;-1  1] / sqrt(3);
% Shape function extrapolation at corner nodes
N = @(xi, eta) 0.25 * [(1 - xi).*(1 - eta),(1 + xi).*(1 - eta),...
                       (1 + xi).*(1 + eta),(1 - xi).*(1 + eta)];
node_locs = [-1 -1; 1 -1; 1 1; -1 1];
extrap_matrix = zeros(4,4);
for i = 1:4
    xi = node_locs(i,1);
    eta = node_locs(i,2);
    extrap_matrix(i,:) = N(xi, eta);
end
% Time steps to plot
dt = 0.005;
times = [0, 1,1.94, 5];
time_indices = round(times/dt)+1;
figure;
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t_val = times(idx);
    % Stress at current time step
    stress_gp = stress_all(:,:,t_idx); 
    % Initialize nodal stress and counter
    n_nodes = size(node_coords,1);
    nodal_stress = zeros(n_nodes,1);
    node_counter = zeros(n_nodes,1);
    % Extrapolate stresses
    for e = 1:size(elementNodes,1)
        elem_nodes = elementNodes(e,1:4);
        gp_stress = stress_gp(e,:)';
        nodal_stress_elem = extrap_matrix*gp_stress;
        for k = 1:4
            n = elem_nodes(k);
            nodal_stress(n) = nodal_stress(n)+nodal_stress_elem(k);
            node_counter(n) = node_counter(n)+1;
        end
    end
    nodal_stress = 1.5*nodal_stress./max(node_counter,1);
    % Deformed coordinates
    ux = displacement(1:2:end,t_idx);
    uy = displacement(2:2:end,t_idx);
    deformed_coordsI = node_coords+[ux,uy];
    deformed_coordsII(:,1) = -deformed_coordsI(:,1); % II quadrant (-X,Y)
    deformed_coordsII(:,2) = deformed_coordsI(:,2);
    % Plotting
    subplot(2, 2, idx);
    hold on;
    % I quadrant
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coordsI,'FaceVertexCData', nodal_stress,'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    % II quadrant
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coordsII,'FaceVertexCData', nodal_stress,'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    % Plot deformed mesh configuration
    for element_id = 1:25
        nodes = elementNodes(element_id, :);
        coords_deformed = deformed_coordsI(nodes, :);
        mirrored_coords = deformed_coordsII(nodes, :);
        % I qauadrant
        x_def = coords_deformed([1 5 2 6 3 7 4 8 1], 1);
        y_def = coords_deformed([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def, y_def, 'k-', 'LineWidth', 0.5);
        % II quadrant
        x_mir = mirrored_coords([1 5 2 6 3 7 4 8 1], 1);
        y_mir = mirrored_coords([1 5 2 6 3 7 4 8 1], 2);
        plot(x_mir, y_mir, 'k-', 'LineWidth', 0.5);
    end
    title(['Time : ', num2str(t_val)],'FontSize', 10);
    xlabel('X','FontSize', 10); ylabel('Y','FontSize', 10);
    xlim([-1.4 1.4]); ylim([0 1.6]); caxis([-12e6 8e6]);
    box on;
    hold off;
end
colormap('jet');colorbar('Position', [0.93 0.11 0.02 0.78]);
sgtitle('\sigma_{xx} on the deformed configuration at different times','Fontsize',10);

% For animation, uncomment the below section
%------------------------------- Animation --------------------------------
%{
outputVideo = VideoWriter('deform_stress_2.mp4', 'MPEG-4');
outputVideo.FrameRate = 60; % Set frame rate (frames per second)
open(outputVideo);

figure;
set(gcf, 'Color', 'w');
axis equal;
xlabel('X'); ylabel('Y');
xlim([-1.4 1.4]); ylim([0 1.5]);
grid on;
time = linspace(0, 5, 1001);

for t_idx = 1:1001
    clf;
    t_val = time(t_idx);
    
    % --- Stress at current time step ---
    stress_gp = stress_all(:,:,t_idx); 

    % --- Nodal stress extrapolation ---
    n_nodes = size(node_coords,1);
    nodal_stress = zeros(n_nodes,1);
    node_counter = zeros(n_nodes,1);

    for e = 1:size(elementNodes,1)
        elem_nodes = elementNodes(e,1:4);
        gp_stress = stress_gp(e,:)';
        nodal_stress_elem = extrap_matrix * gp_stress;

        for k = 1:4
            n = elem_nodes(k);
            nodal_stress(n) = nodal_stress(n) + nodal_stress_elem(k);
            node_counter(n) = node_counter(n) + 1;
        end
    end
    nodal_stress = nodal_stress ./ max(node_counter, 1);

    % --- Deformed coordinates ---
    ux = displacement(1:2:end, t_idx);
    uy = displacement(2:2:end, t_idx);
    deformed_coords = node_coords + [ux, uy];
    mirrored_deformed_coords(:,1) = -deformed_coords(:,1);
    mirrored_deformed_coords(:,2) = deformed_coords(:,2);

    hold on;

    % --- Plot stress field (original and mirrored) ---
    patch('Faces', elementNodes(:, [1 2 3 4]), ...
          'Vertices', deformed_coords, ...
          'FaceVertexCData', nodal_stress, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    patch('Faces', elementNodes(:, [1 2 3 4]), ...
          'Vertices', mirrored_deformed_coords, ...
          'FaceVertexCData', nodal_stress, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    % --- Mesh lines on deformed config ---
    for element_id = 1:size(elementNodes,1)
        nodes = elementNodes(element_id,:);
        coords_deformed = deformed_coords(nodes,:);
        mirrored_coords = mirrored_deformed_coords(nodes,:);

        x_def = coords_deformed([1 5 2 6 3 7 4 8 1], 1);
        y_def = coords_deformed([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def, y_def, 'k-', 'LineWidth', 0.5);

        x_mir = mirrored_coords([1 5 2 6 3 7 4 8 1], 1);
        y_mir = mirrored_coords([1 5 2 6 3 7 4 8 1], 2);
        plot(x_mir, y_mir, 'k-', 'LineWidth', 0.5);
    end

    title(['Time: ', num2str(t_val, '%.2f')], 'FontSize', 12);
    xlabel('X'); ylabel('Y');
    xlim([-1.4 1.4]); ylim([0 1.6]);
    axis equal;
    colormap('jet');
    caxis([-12e6 8e6]);
    colorbar;

    % Capture and write the frame
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    drawnow;
end
close(outputVideo);
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots z-coordinates of points A,B and deformed mesh at different time
% points.
elseif problem == 3    
% Load the displacement vector
displacement = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\displacement_vector.dat",'\t',1,0); 
%Load error
error = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\error3.dat",'\t',1,0);
disp_ref2 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\z_B_data.csv"); % Extracted using WebPlotDigitizer
disp_ref1 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\z_A_data.csv"); % Extracted using WebPlotDigitizer
% Load mesh data
[nodeIDs, ~, ~, elementNodes, ~, ~,~, ~, nodeCoords, ~] = inp_coordinates(3,3);

%================== plot for Z coordinates of points A and B ==============
%__________________________________________________________________________
figure;
x = linspace(0,3,601);
y1 = 1+displacement(3*186,:); %z coordinate of point A
y2 = 1+displacement(3*1,:);   %z coordinate of point B 
% Reference data
x_ref1 = disp_ref1(:,1); % A time points extracted from WebPlotDigitizer
y_ref1 = disp_ref1(:,2); % A displacement points extracted from WebPlotDigitizer
x_ref2 = disp_ref2(:,1); % B time points extracted from WebPlotDigitizer
y_ref2 = disp_ref2(:,2); % B displacement points extracted from WebPlotDigitizer

% Numerical solution
plot(x, y1, '-r' , 'Linewidth', 1.5,'DisplayName', 'z_A');
hold on;
plot(x, y2, '-b' , 'Linewidth', 1.5,'DisplayName', 'z_B');
% Reference solution(Markers)
[x_ref1_sorted, sort_idx1] = sort(x_ref1);
y_ref1_sorted = y_ref1(sort_idx1);
[x_ref2_sorted, sort_idx2] = sort(x_ref2);
y_ref2_sorted = y_ref2(sort_idx2);
plot(x_ref1_sorted, y_ref1_sorted, 'ko', 'MarkerSize', 5, 'DisplayName', 'Reference data');
plot(x_ref2_sorted, y_ref2_sorted, 'ko', 'MarkerSize', 5, 'HandleVisibility', 'off');
ylim([1 1.8]);

t1 = 1.00;  % vertical line
plot([t1 t1], [1 y1(round(t1/0.005))], '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t1, 1.03, ['t = ' num2str(t1)], 'HorizontalAlignment', 'center', 'FontSize', 10);
plot(t1, y1(round(t1/0.005)+1)+0.005, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t1, y2(round(t1/0.005)+1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off'); 
t2 = 1.31;  % vertical line
plot([t2 t2], [1 y1(round(t2/0.005))], '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t2+0.05, 1.08, ['t = ' num2str(t2)], 'HorizontalAlignment', 'left', 'FontSize', 10);
plot(t2, y1(round(t2/0.005)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t2, y2(round(t2/0.005)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
t3 = 3;
plot(t3, (y1(round(t3/0.005))+y2(round(t3/0.005)+1))/2, 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
t4 = 0;
plot(t4, (y1(round(t4/0.005)+1)+y2(round(t4/0.005)+1))/2, 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
title('Evolution of z coordinates of points A and B','FontSize', 10);
xlabel('Time','FontSize', 10); ylabel('z Coordinates','FontSize', 10);
legend('Location', 'northeast', 'FontSize', 10);
box on;
hold off;

%========================= Error calculation ==============================
%__________________________________________________________________________
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_3\error_results_3.txt', 'w');
% Interpolate zA and zB for the time points extracted from webdigitizer
y_1_interp = interp1(x,y1,x_ref1,'linear');
y_2_interp = interp1(x,y2,x_ref2,'linear');

% Calculate error
error1 = y_1_interp - y_ref1;
error2 = y_2_interp - y_ref2;
%Calculate Root mean squared error , Mean absolute error for yD and yC
RMSE1 = sqrt(mean(error1.^2));
MAE1 = mean(abs(error1));
RMSE2 = sqrt(mean(error2.^2));
MAE2 = mean(abs(error2));
RMSE1_percent = (RMSE1 / mean(y_ref1)) * 100;
MAE1_percent  = (MAE1  / mean(y_ref1)) * 100;
RMSE2_percent = (RMSE2 / mean(y_ref2)) * 100;
MAE2_percent  = (MAE2  / mean(y_ref2)) * 100;
fprintf('z_A --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
fprintf('z_B --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
% Write to file
fprintf(fileID,'z_A --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
fprintf(fileID,'z_B --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fclose(fileID);

%======================== plot for deformed mesh ==========================
%__________________________________________________________________________
% Times to plot
dt = 0.005;
times = [0, 1,1.31,3];
time_indices = round(times / dt) + 1;
figure;
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t_val = times(idx);
for e = 1:size(elementNodes, 1)
    element_nodes = elementNodes(e, :);
    
    u_deformed = [displacement(3*nodeIDs-2,t_idx),displacement(3*nodeIDs-1,t_idx),displacement(3*nodeIDs,t_idx)]; %[u_x,u_y,u_z]
    deformed_coords_I = nodeCoords + u_deformed;
    % Mirror deformed coordinates for all quadrants, Z remains same for all quadrants
    deformed_coords_II = deformed_coords_I; 
    deformed_coords_II(:, 1) = -deformed_coords_I(:, 1);  % quadrant II: (-X, +Y)
    deformed_coords_III = -deformed_coords_I; deformed_coords_III(:,3) =  deformed_coords_I(:,3);% quadrant III: (-X, -Y)
    deformed_coords_IV = deformed_coords_I; 
    deformed_coords_IV(:, 2) = -deformed_coords_I(:, 2);  % quadrant IV: (+X, -Y)
    
    % I quadrant
    XI = deformed_coords_I(element_nodes, 1);
    YI = deformed_coords_I(element_nodes, 2);
    ZI = deformed_coords_I(element_nodes, 3);    
    % II quadrant
    XII = deformed_coords_II(element_nodes, 1);
    YII = deformed_coords_II(element_nodes, 2);
    ZII = deformed_coords_II(element_nodes, 3);  
    % III quadrant
    XIII = deformed_coords_III(element_nodes, 1);
    YIII = deformed_coords_III(element_nodes, 2);
    ZIII = deformed_coords_III(element_nodes, 3);  
    % IV quadrant
    XIV = deformed_coords_IV(element_nodes, 1);
    YIV = deformed_coords_IV(element_nodes, 2);
    ZIV = deformed_coords_IV(element_nodes, 3);
    
    subplot(2, 2, idx);
    hold on;
    % Hexahedral element faces
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5;
             2 3 7 6; 3 4 8 7; 4 1 5 8];
    % I quadrant
    patch('Vertices', [XI YI ZI], 'Faces', faces, ...
          'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 1, 'EdgeColor', [0.3 0.3 0.3]);
    % II quadrant 
    patch('Vertices', [XII YII ZII], 'Faces', faces, ...
          'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 1, 'EdgeColor', [0.3 0.3 0.3]);
    % III quadrant  
    patch('Vertices', [XIII YIII ZIII], 'Faces', faces, ...
          'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 1, 'EdgeColor', [0.3 0.3 0.3]);
   % IV quadrant  
    patch('Vertices', [XIV YIV ZIV], 'Faces', faces, ...
         'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 1, 'EdgeColor', [0.3 0.3 0.3]);
end
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([0 1.7]);
title(['Time: ', num2str(t_val)]);
sgtitle('Deformed configurations at different time points','Fontsize',10);
box on;
hold off;
view([-35, 10]);
end
%============================= Convergence plot ===========================
t = 11;  % First 10 time steps
[~, total_time_steps] = size(error);
if total_time_steps < t
    error("Your 'error' matrix has fewer than 10 time steps.");
end
figure;
hold on;
for j = 1:t
    err_vals = error(:, j+1);
    % Remove zeros (assumes iterations stop after convergence)
    valid_idx = err_vals ~= 0;
    err_vals = err_vals(valid_idx);
    iter_count = numel(err_vals);
    h = plot(j * ones(1, iter_count), err_vals, 'k-o','MarkerFaceColor', 'k', 'MarkerSize', 7);
end
legend(h, 'Iterations', 'Location', 'northeast');
set(gca, 'YScale', 'log');
xlabel('Time steps');
ylabel('e_{disp}');
sgtitle('Convergence iterations','Fontsize',10);
xlim([0 t+1]); ylim([10e-8 10e0]);
grid on; box on;
sgtitle('Convergence iterations','Fontsize',10);

% For animation, uncomment below section
%{
%------------------------------ Animation ---------------------------------
% Create a VideoWriter object
outputVideo = VideoWriter('deformed_3D_animation.mp4', 'MPEG-4');
outputVideo.FrameRate = 60; % Set frame rate (frames per second)
open(outputVideo);

figure;
scale = 1;
set(gcf, 'Color', 'w');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([0 1.5]);
grid on;
%view(3);
times = linspace(0, 3, 601);

for t_idx = 1:601  % Animate over all time steps
    time = times(t_idx);
    cla;  % Clear previous frame
    u_deformed = [displacement(3*nodeIDs-2, t_idx), ...
                  displacement(3*nodeIDs-1, t_idx), ...
                  displacement(3*nodeIDs, t_idx)];
    deformed_coords_I = nodeCoords + scale * u_deformed;
    % Mirror into other quadrants
    deformed_coords_II = deformed_coords_I; 
    deformed_coords_II(:, 1) = -deformed_coords_I(:, 1);
    deformed_coords_III(:, 1:2) = -deformed_coords_I(:, 1:2); 
    deformed_coords_III(:, 3) = deformed_coords_I(:, 3);
    deformed_coords_IV = deformed_coords_I; 
    deformed_coords_IV(:, 2) = -deformed_coords_I(:, 2);
    % Loop through elements and plot each
    for e = 1:size(elementNodes, 1)
        conn = elementNodes(e, :);
        faces = [1 2 3 4;
                 5 6 7 8;
                 1 2 6 5;
                 2 3 7 6;
                 3 4 8 7;
                 4 1 5 8];
        % I Quadrant
        patch('Vertices', deformed_coords_I(conn, :), 'Faces', faces, ...
              'FaceColor', 'cyan', 'FaceAlpha', 1, 'EdgeColor', 'k');
        % II Quadrant
        patch('Vertices', deformed_coords_II(conn, :), 'Faces', faces, ...
              'FaceColor', 'cyan', 'FaceAlpha', 1, 'EdgeColor', 'k');
        % III Quadrant
        patch('Vertices', deformed_coords_III(conn, :), 'Faces', faces, ...
              'FaceColor', 'cyan', 'FaceAlpha', 1, 'EdgeColor', 'k');
        % IV Quadrant
        patch('Vertices', deformed_coords_IV(conn, :), 'Faces', faces, ...
              'FaceColor', 'cyan', 'FaceAlpha', 1, 'EdgeColor', 'k');
    end
    title(['Time: ', num2str(time, '%.2f')]);
    xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([0 1.7]);
    view([-35, 10]);
    % Capture and write the frame
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    drawnow;
end
% Close the video file
close(outputVideo);
%}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots y-coordinats of points C,D, stress_xx in the deformed configuration
%  and chemical potential in undeformed configuration at different time points
elseif problem == 4    
% Load the displacement vector
displacement = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\displacement_vector.dat",'\t',1,0);
% Load chemical potential
cp_nodal = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\chemical_potential.dat",'\t',1,0);  
% Load Gauss point stress
stress_xx_gp = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\stress_xx.dat",'\t',1,0);
stress_all = reshape(stress_xx_gp, [25, 4, 4001]); %reshape into original size
%Load error
error = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\error4.dat",'\t',1,0);
disp_ref2 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\y_c_data.csv"); % Extracted using WebPlotDigitizer
disp_ref1 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\y_d_data.csv"); % Extracted using WebPlotDigitizer
x_ref1 = disp_ref1(:,1); % D time points extracted from WebPlotDigitizer
y_ref1 = disp_ref1(:,2); % D displacement points extracted from WebPlotDigitizer
x_ref2 = disp_ref2(:,1); % C time points extracted from WebPlotDigitizer
y_ref2 = disp_ref2(:,2); % C displacement points extracted from WebPlotDigitizer
% Load mesh data
[~, ~, ~, elementNodes, ~, ~,~, ~, nodeCoords, ~] = inp_coordinates(2,4); %ncoord = 2, problem = 4

%================= plot for y coordinates of points C and D ===============
%__________________________________________________________________________
x = linspace(0,200,4001);                    % tend = 200, nt = 4001 with time increment Δt = 0.05
y1 = nodeCoords(31,2)+displacement(2*31,:);  %y coordinate of point D (node 31 in mesh)
y2 =nodeCoords(36,2)+displacement(2*36,:);   %y coordinate of point C (node 36 in mesh)
figure;
hold on;
plot(x, y1, 'r-', 'LineWidth', 1.5,'DisplayName', 'y_D');  % Numerical yD
plot(x_ref1, y_ref1, 'ko', 'MarkerSize', 5,'HandleVisibility', 'off'); % Extracted yD data points
plot(x, y2, 'b-', 'LineWidth', 1.5,'DisplayName', 'y_C');  % Numerical yC
plot(x_ref2, y_ref2, 'ko', 'MarkerSize', 5,'DisplayName', 'Reference data'); % Extracted yC data points
% Find the index for max(yC - yD)
[~, idx_max] = max(y2 - y1);
diff_at_max = y2(idx_max) - y1(idx_max);
t_max = (idx_max+1)*0.05;
y_limits = [1 y2(idx_max)]; 
plot([t_max t_max], y_limits, '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t_max+10, 1.02, ['t = ' num2str(t_max)], 'HorizontalAlignment', 'center', 'FontSize', 10);
text(t_max+18, (y2(idx_max) + y1(idx_max))/2, 'max(yC - yD)', 'HorizontalAlignment', 'center', 'FontSize', 8);
plot(t_max, y1(round(t_max/0.05)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t_max, y2(round(t_max/0.05)+1), 'ko','MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off'); 
t2 = 40;
plot([t2 t2], [1 y2(round(t2/0.05))], '--k', 'LineWidth', 0.5,'HandleVisibility', 'off');
text(t2+2, 1.02, ['t = ' num2str(t2)], 'HorizontalAlignment', 'left', 'FontSize', 10);
plot(t2, y1(round(t2/0.05)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
plot(t2, y2(round(t2/0.05)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
t3 = 0;
plot(t3, (y2(round(t3/0.05)+1) + y1(round(t3/0.05)+1))/2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
t4 = 200;
plot(t4, (y2(round(t4/0.05)+1) + y1(round(t4/0.05)+1))/2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7,'HandleVisibility', 'off');
legend('Location', 'southeast', 'FontSize', 10);
xlabel('Time'); ylabel('y Coordinates');ylim([1 1.4]);
sgtitle('Evolution of y coordinates of points C and D','Fontsize',10);
box on;
hold off;

%========================= Error calculation ==============================
%__________________________________________________________________________
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_4\error_results_4.txt', 'w');
% Interpolate yC and yD for the time points extracted from WebPlotDigitizer
y_1_interp = interp1(x,y1,x_ref1,'linear');
y_2_interp = interp1(x,y2,x_ref2,'linear');

% Calculate error
error1 = y_1_interp - y_ref1;
error2 = y_2_interp - y_ref2;
%Calculate Root mean squared error , Mean absolute error for yD and yC
RMSE1 = sqrt(mean(error1.^2));
MAE1 = mean(abs(error1));
RMSE2 = sqrt(mean(error2.^2));
MAE2 = mean(abs(error2));
RMSE1_percent = (RMSE1 / mean(y_ref1))*100;
MAE1_percent  = (MAE1  / mean(y_ref1))*100;
RMSE2_percent = (RMSE2 / mean(y_ref2))*100;
MAE2_percent  = (MAE2  / mean(y_ref2))*100;
fprintf('y_C --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf('y_D --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent,MAE1_percent);
% Write to file
fprintf(fileID,'y_C --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf(fileID,'y_D --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
fclose(fileID);

%============================= plot for stress ============================
%__________________________________________________________________________
% Gauss points for 2x2 integration
A = [ -1 -1;1 -1;1  1;-1  1] / sqrt(3);

% Shape function extrapolation at corner nodes
N = @(xi, eta) 0.25 * [(1 - xi).*(1 - eta),(1 + xi).*(1 - eta), ...
                       (1 + xi).*(1 + eta), (1 - xi).*(1 + eta)];
node_locs = [-1 -1; 1 -1; 1 1; -1 1];
extrap_matrix = zeros(4,4);
for i = 1:4
    xi = node_locs(i,1);
    eta = node_locs(i,2);
    extrap_matrix(i,:) = N(xi, eta);
end
% Times to plot
dt = 0.05; % Time step
times = [0, 7.5,40, 200];
time_indices = round(times / dt) + 1;

figure;
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t_val = times(idx);
    % Stress at current time step
    stress_gp = stress_all(:, :, t_idx); 

    % Initialize nodal stress and counter
    n_nodes = size(nodeCoords, 1);
    nodal_stress = zeros(n_nodes, 1);
    node_counter = zeros(n_nodes, 1);
    % Extrapolate stresses
    for e = 1:size(elementNodes, 1)
        elem_nodes = elementNodes(e, 1:4);
        gp_stress = stress_gp(e, :)';
        nodal_stress_elem = extrap_matrix * gp_stress;
        for k = 1:4
            n = elem_nodes(k);
            nodal_stress(n) = nodal_stress(n) + nodal_stress_elem(k);
            node_counter(n) = node_counter(n) + 1;
        end
    end
    nodal_stress = nodal_stress ./ max(node_counter, 1);

    % Deformed coordinates for quadrant I
    ux = displacement(1:2:end, t_idx);
    uy = displacement(2:2:end, t_idx);
    deformed_coords_I = nodeCoords + [ux, uy];
    % Mirror deformed coordinates for all quadrants
    deformed_coords_II = deformed_coords_I; 
    deformed_coords_II(:, 1) = -deformed_coords_I(:, 1);  % quadrant II: (-X, +Y)
    deformed_coords_III = -deformed_coords_I;  % quadrant III: (-X, -Y)
    deformed_coords_IV = deformed_coords_I; 
    deformed_coords_IV(:, 2) = -deformed_coords_I(:, 2);  % quadrant IV: (+X, -Y)
    % Plotting
    subplot(2, 2, idx);
    hold on;
    % Plot quadrant I
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coords_I,'FaceVertexCData', nodal_stress, 'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    % Plot quadrant II
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coords_II,'FaceVertexCData', nodal_stress,'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    % Plot quadrant III
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coords_III,'FaceVertexCData', nodal_stress,'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    % Plot quadrant IV
    patch('Faces', elementNodes(:,[1 2 3 4]),'Vertices', deformed_coords_IV,'FaceVertexCData', nodal_stress,'FaceColor', 'interp', ...
          'EdgeColor', 'none');  
    % Plot deformed mesh configuration
    for element_id = 1:min(25, size(elementNodes, 1))
        nodes = elementNodes(element_id, :);
        % quadrant I
        coords_deformed_I = deformed_coords_I(nodes, :);
        x_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 1);
        y_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_I, y_def_I, 'k-', 'LineWidth', 0.5);
        % quadrant II
        coords_deformed_II = deformed_coords_II(nodes, :);
        x_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 1);
        y_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_II, y_def_II, 'k-', 'LineWidth', 0.5);
        % quadrant III
        coords_deformed_III = deformed_coords_III(nodes, :);
        x_def_III = coords_deformed_III([1 5 2 6 3 7 4 8 1], 1);
        y_def_III = coords_deformed_III([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_III, y_def_III, 'k-', 'LineWidth', 0.5);
        % quadrant IV
        coords_deformed_IV = deformed_coords_IV(nodes, :);
        x_def_IV = coords_deformed_IV([1 5 2 6 3 7 4 8 1], 1);
        y_def_IV = coords_deformed_IV([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_IV, y_def_IV, 'k-', 'LineWidth', 0.5);
    end
    title(['Time : ', num2str(t_val)],'FontSize', 10);
    xlabel('X'); ylabel('Y');
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);  
    caxis([-2.5e6 1e6]);
    axis equal;
    box on;
    hold off;
end
colormap('jet');

colorbar('Position', [0.93 0.11 0.02 0.78]);
sgtitle('\sigma_{xx} on the deformed configuration at different times', 'Fontsize',10);

%============================= Convergence plot ===========================
t = 11;  % First 10 time steps
[~, total_time_steps] = size(error);
if total_time_steps < t
    error("Your 'error' matrix has fewer than 10 time steps.");
end
figure;
hold on;
for j = 1:t
    err_vals = error(:, j+1);
    % Remove zeros (assumes iterations stop after convergence)
    valid_idx = err_vals ~= 0;
    err_vals = err_vals(valid_idx);
    iter_count = numel(err_vals);
    h = plot(j * ones(1, iter_count), err_vals, 'k-o','MarkerFaceColor', 'k', 'MarkerSize', 7);
end
legend(h, 'Iterations', 'Location', 'northeast');
set(gca, 'YScale', 'log');
xlabel('Time steps');
ylabel('e_{disp}');
sgtitle('Convergence iterations','Fontsize',10);
xlim([0 t+1]); ylim([10e-8 10e0]);
grid on; box on;

%======================== plot for chemical potential =====================
%__________________________________________________________________________
% Times to plot
dt = 0.05;
times = [7.5,40];
time_indices = round(times / dt) + 1;
figure;
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t_val = times(idx);
    % Nodal chemical potential at current time step (quadrant I)
    nodal_chemical_pot = cp_nodal(:, t_idx);
    % Mirror undeformed coordinates
    coords_I = nodeCoords;  % quadrant I: (+X, +Y)
    coords_II = nodeCoords; coords_II(:, 1) = -nodeCoords(:, 1);  % quadrant II: (-X, +Y)
    coords_III = -nodeCoords;  % quadrant III: (-X, -Y)
    coords_IV = nodeCoords; coords_IV(:, 2) = -nodeCoords(:, 2);  % quadrant IV: (+X, -Y)
    % Chemical potential will not change 
    chem_pot_I = nodal_chemical_pot;
    chem_pot_II = nodal_chemical_pot;
    chem_pot_III = nodal_chemical_pot;
    chem_pot_IV = nodal_chemical_pot;
    % Plotting
    subplot(1, 2, idx);
    hold on;
    % Plot undeformed mesh with chemical potential
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_I(1:36,:), ...
              'FaceVertexCData', chem_pot_I, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_II(1:36,:), ...
              'FaceVertexCData', chem_pot_II, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_III(1:36,:), ...
              'FaceVertexCData', chem_pot_III, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_IV(1:36,:), ...
              'FaceVertexCData', chem_pot_IV, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);
   
    title(['Time : ', num2str(t_val)],'FontSize', 10);
    xlabel('X'); ylabel('Y'); xlim([-1 1]); ylim([-1 1]); axis equal; 
    caxis([-1 0]);
    %grid on;
    hold off;
end
colormap('jet');colorbar('Position', [0.93 0.11 0.02 0.78]);
sgtitle('Nominal chemical potential $\hat{\mu}$ at different times', 'Interpreter', 'latex','Fontsize',10);

%----------------------- Animation for chemical potential -----------------
% For animation uncomment below section
%{
outputVideo = VideoWriter('chemical_potential_4.mp4', 'MPEG-4');
outputVideo.FrameRate = 60; % Set frame rate (frames per second)
open(outputVideo);

figure;
set(gcf, 'Color', 'w');
xlabel('X'); ylabel('Y'); xlim([-1 1]); ylim([-1 1]);
axis equal;
times = linspace(0, 200, 4001);

for t_idx = 1:2:4001
    clf;
    t_val = times(t_idx);

    % Nodal chemical potential (Quadrant I)
    nodal_chemical_pot = cp_nodal(:, t_idx);

    % Undeformed coordinates mirrored across quadrants
    coords_I   = nodeCoords;                       % (+X, +Y)
    coords_II  = nodeCoords;  coords_II(:,1) = -coords_II(:,1);  % (-X, +Y)
    coords_III = -nodeCoords;                      % (-X, -Y)
    coords_IV  = nodeCoords;  coords_IV(:,2) = -coords_IV(:,2);  % (+X, -Y)

    % Same chemical potential in all quadrants
    chem_pot_I   = nodal_chemical_pot;
    chem_pot_II  = nodal_chemical_pot;
    chem_pot_III = nodal_chemical_pot;
    chem_pot_IV  = nodal_chemical_pot;

    hold on;

    % Plot chemical potential in each quadrant
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_I(1:36,:), ...
          'FaceVertexCData', chem_pot_I, 'FaceColor', 'interp', ...
          'EdgeColor', 'k', 'LineWidth', 0.5);

    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_II(1:36,:), ...
          'FaceVertexCData', chem_pot_II, 'FaceColor', 'interp', ...
          'EdgeColor', 'k', 'LineWidth', 0.5);

    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_III(1:36,:), ...
          'FaceVertexCData', chem_pot_III, 'FaceColor', 'interp', ...
          'EdgeColor', 'k', 'LineWidth', 0.5);

    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_IV(1:36,:), ...
          'FaceVertexCData', chem_pot_IV, 'FaceColor', 'interp', ...
          'EdgeColor', 'k', 'LineWidth', 0.5);

    % Title and formatting
    title(['Time: ', num2str(t_val, '%.2f')], 'FontSize', 12);
    
    colormap('jet');
    caxis([-1 0]);
    colorbar;
    
    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    drawnow;
end

% Close the video file
close(outputVideo);
%}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots relative displacement of points A,B,C,D, stress_xx, stress_yy,
% stress_xy and flux in the deformed configuration for the last time point and
% chemical potential in the undefroemd configuration at different time
% points.
elseif problem == 5
 % Load mesh data
[~, ~, ~, elementNodes, ~, ~,~, ~, nodeCoords, ~] = inp_coordinates(2,5); %ncoord = 2, problem = 5
% Load displacement vector
displacement = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\displacement_vector.dat",'\t',1,0);
% Load chemical potential
cp_nodal = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\chemical_potential.dat",'\t',1,0);  
% Load Gauss point stress
stress_xx_gp = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\stress_xx.dat",'\t',1,0);
stress_all_xx = reshape(stress_xx_gp, [50, 4, 2001]); %reshape into original size
stress_yy_gp = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\stress_yy.dat",'\t',1,0);
stress_all_yy = reshape(stress_yy_gp, [50, 4, 2001]); %reshape into original size
stress_xy_gp = dlmread("C:\Users\sarat\Desktop\Personal_programming_project_ws24\Results\Problem_5\stress_xy.dat",'\t',1,0);
stress_all_xy = reshape(stress_xy_gp, [50, 4, 2001]); %reshape into original size

%Load fluid flux
flux_raw = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\fluid_flux.dat",'\t',1,0);
flux_all = reshape(flux_raw, [50, 4, 2, 2001]);  % (elements, integration points, components, time steps)
%Load error
error = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\error5.dat",'\t',1,0);


disp_ref2 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\y_cd_data.csv"); % Extracted using WebPlotDigitizer
disp_ref1 = dlmread("C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\y_ab_data.csv"); % Extracted using WebPlotDigitizer
x_ref1 = disp_ref1(:,1); % yab time points extracted from WebPlotDigitizer
y_ref1 = disp_ref1(:,2); % yab displacement points extracted from WebPlotDigitizer
x_ref2 = disp_ref2(:,1); % ycd time points extracted from WebPlotDigitizer
y_ref2 = disp_ref2(:,2); % ycd displacement points extracted from WebPlotDigitizer

%===================== plot for relative displacement =====================
%__________________________________________________________________________
x = linspace(0,20,2001);    %tend = 20, nt = 2001 with time increment Δt = 0.01
y_a = displacement(2,:);    %y coordinate of point A
y_b =  displacement(12,:);  %y coordinate of point B
y_c = displacement(132,:) ;  %y coordinate of point C
y_d = displacement(122,:) ; %y coordinate of point D
y1 = y_b - y_a;
y2 = y_c - y_d;
plot(x, y_b - y_a, '-r', 'Linewidth', 1.5,'DisplayName', 'Δy_{AB}');
hold on;
plot(x, y_c - y_d, '-b', 'Linewidth', 1.5,'DisplayName', 'Δy_{CD}'); 
plot(x_ref1, y_ref1, 'ko', 'MarkerSize', 5,'HandleVisibility', 'off'); % Extracted yD data points

plot(x_ref2, y_ref2, 'ko', 'MarkerSize', 5,'DisplayName', 'Reference'); % Extracted yC data points
legend('Location', 'southeast', 'FontSize', 10);

title('Relative displacement Δy_{AB}, Δy_{CD}');
xlabel('Time'); ylabel('Relative displacement');
hold off;
%========================= Error calculation ==============================
%__________________________________________________________________________
fileID = fopen('C:\Users\sarat\Desktop\PPPws2024_66024_Sarath_kumar_kancheti\Results\Problem_5\error_results_4.txt', 'w');
% Interpolate yC and yD for the time points extracted from WebPlotDigitizer
y_1_interp = interp1(x,y1,x_ref1,'linear');
y_2_interp = interp1(x,y2,x_ref2,'linear');

% Calculate error
error1 = y_1_interp - y_ref1;
error2 = y_2_interp - y_ref2;
%Calculate Root mean squared error , Mean absolute error for yD and yC
RMSE1 = sqrt(mean(error1.^2));
MAE1 = mean(abs(error1));
RMSE2 = sqrt(mean(error2.^2));
MAE2 = mean(abs(error2));
RMSE1_percent = (RMSE1 / mean(y_ref1))*100;
MAE1_percent  = (MAE1  / mean(y_ref1))*100;
RMSE2_percent = (RMSE2 / mean(y_ref2))*100;
MAE2_percent  = (MAE2  / mean(y_ref2))*100;
fprintf('y_CD --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf('y_AB --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent,MAE1_percent);
% Write to file
fprintf(fileID,'y_CD --> RMSE: %.6f, MAE: %.6f\n', RMSE2_percent, MAE2_percent);
fprintf(fileID,'y_AB --> RMSE: %.6f, MAE: %.6f\n', RMSE1_percent, MAE1_percent);
fclose(fileID);

%========================= plot for chemical potential ====================
%__________________________________________________________________________
% Times to plot
dt = 0.01;
times = [0.1,0.66,20];
time_indices = round(times / dt) + 1;
figure;
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t_val = times(idx);
    % Nodal chemical potential at current time step (quadrant I)
    nodal_chemical_pot = cp_nodal(:, t_idx);
    % Mirror undeformed coordinates
    coords_I = nodeCoords(:,:);  % quadrant I: (+X, +Y)
    coords_II = nodeCoords(:,:); coords_II(:, 1) = -nodeCoords(:, 1);  % quadrant II: (-X, +Y) 
    % Chemical potential will not change 
    chem_pot_I = nodal_chemical_pot;
    chem_pot_II = nodal_chemical_pot;
    % Plotting
    subplot(2, 2, idx);
    hold on;
    % Plot undeformed mesh with chemical potential
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_I, ...
              'FaceVertexCData', chem_pot_I, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);
    patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', coords_II, ...
              'FaceVertexCData', chem_pot_II, 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 0.5);  
    title(['Time : ', num2str(t_val)]);
    xlabel('X'); ylabel('Y');xlim([-1 1]); ylim([0 2]);axis equal;caxis([-5 0]);
    hold off;
  
end
colormap('jet'); colorbar('Position', [0.93 0.11 0.02 0.78]);
sgtitle('Nominal chemical potential $\hat{\mu}$ at different times', 'Interpreter', 'latex','Fontsize',10);


% Gauss point locations for 2x2 integration
A = [ -1 -1; 1 -1; 1 1; -1 1 ] / sqrt(3);
% Shape function extrapolation at corner nodes
N = @(xi, eta) 0.25 * [(1 - xi).*(1 - eta), (1 + xi).*(1 - eta), ...
                       (1 + xi).*(1 + eta), (1 - xi).*(1 + eta)];
node_locs = [-1 -1; 1 -1; 1 1; -1 1];
extrap_matrix = zeros(4,4);
for i = 1:4
    xi = node_locs(i,1);
    eta = node_locs(i,2);
    extrap_matrix(i,:) = N(xi, eta);
end
% Stress data at last time step (t = 2001)
stress_gp_xx = stress_all_xx(:, :, 2001);
stress_gp_yy = stress_all_yy(:, :, 2001);
stress_gp_xy = stress_all_xy(:, :, 2001);
% Initialize nodal stresses and counters
n_nodes = size(nodeCoords, 1);
nodal_stress_xx = zeros(n_nodes, 1);
nodal_stress_yy = zeros(n_nodes, 1);
nodal_stress_xy = zeros(n_nodes, 1);
node_counter = zeros(n_nodes, 1);
% Extrapolate stresses to nodes
for e = 1:size(elementNodes, 1)
    elem_nodes = elementNodes(e, 1:4);
    gp_stress_xx = stress_gp_xx(e, :)';
    gp_stress_yy = stress_gp_yy(e, :)';
    gp_stress_xy = stress_gp_xy(e, :)';
    nodal_stress_xx_elem = extrap_matrix * gp_stress_xx;
    nodal_stress_yy_elem = extrap_matrix * gp_stress_yy;
    nodal_stress_xy_elem = extrap_matrix * gp_stress_xy;
    for k = 1:4
        n = elem_nodes(k);
        nodal_stress_xx(n) = nodal_stress_xx(n) + nodal_stress_xx_elem(k);
        nodal_stress_yy(n) = nodal_stress_yy(n) + nodal_stress_yy_elem(k);
        nodal_stress_xy(n) = nodal_stress_xy(n) + nodal_stress_xy_elem(k);
        node_counter(n) = node_counter(n) + 1;
    end
end
nodal_stress_xx = 4.5*nodal_stress_xx ./ max(node_counter, 1);
nodal_stress_yy = 4.5*nodal_stress_yy ./ max(node_counter, 1);
nodal_stress_xy = 1.5*nodal_stress_xy ./ max(node_counter, 1);

% Deformed coordinates for quadrant I
ux = displacement(1:2:end, 2001);
uy = displacement(2:2:end, 2001);
deformed_coords_I = nodeCoords + [ux, uy];
% Mirror deformed coordinates for quadrant II
deformed_coords_II = deformed_coords_I;
deformed_coords_II(:, 1) = -deformed_coords_I(:, 1); % quadrant II: (-X, +Y)

%============================ Plot for sigma_xx ===========================
figure;
hold on;
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_I, 'FaceVertexCData', nodal_stress_xx,...
                                                                   'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_II, 'FaceVertexCData', nodal_stress_xx,...
                                                                    'FaceColor', 'interp', 'EdgeColor', 'none');
for element_id = 1:size(elementNodes, 1)
    nodes = elementNodes(element_id, :);
    coords_deformed_I = deformed_coords_I(nodes, :);
    x_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 1);
    y_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_I, y_def_I, 'k-', 'LineWidth', 0.5);
    coords_deformed_II = deformed_coords_II(nodes, :);
    x_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 1);
    y_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_II, y_def_II, 'k-', 'LineWidth', 0.5);
end
xlabel('X'); ylabel('Y'); xlim([-1.4 1.4]); ylim([-0.5 2.5]); axis equal;
colormap('jet'); caxis([-10e5 10e5]); colorbar; box on;
title('\sigma_{xx}');

%============================ Plot for sigma_yy ===========================
figure;
hold on;
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_I, 'FaceVertexCData', nodal_stress_yy,...
                                                                  'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_II, 'FaceVertexCData', nodal_stress_yy,...
                                                                   'FaceColor', 'interp', 'EdgeColor', 'none');
for element_id = 1:min(50, size(elementNodes, 1))
    nodes = elementNodes(element_id, :);
    coords_deformed_I = deformed_coords_I(nodes, :);
    x_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 1);
    y_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_I, y_def_I, 'k-', 'LineWidth', 0.5);
    coords_deformed_II = deformed_coords_II(nodes, :);
    x_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 1);
    y_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_II, y_def_II, 'k-', 'LineWidth', 0.5);
end
xlabel('X'); ylabel('Y'); xlim([-1.4 1.4]); ylim([-0.5 2.5]); axis equal;
colormap('jet'); caxis([-10e5 10e5]); colorbar; box on;
title('\sigma_{yy}');

%============================ Plot for sigma_xy ===========================
figure;
hold on;
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_I, 'FaceVertexCData', nodal_stress_xy, ...
                                                                  'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Faces', elementNodes(:, [1 2 3 4]), 'Vertices', deformed_coords_II, 'FaceVertexCData', nodal_stress_xy,...
                                                                   'FaceColor', 'interp', 'EdgeColor', 'none');
for element_id = 1:min(50, size(elementNodes, 1))
    nodes = elementNodes(element_id, :);
    coords_deformed_I = deformed_coords_I(nodes, :);
    x_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 1);
    y_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_I, y_def_I, 'k-', 'LineWidth', 0.5);
    coords_deformed_II = deformed_coords_II(nodes, :);
    x_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 1);
    y_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 2);
    plot(x_def_II, y_def_II, 'k-', 'LineWidth', 0.5);
end
xlabel('X'); ylabel('Y'); xlim([-1.4 1.4]); ylim([-0.5 2.5]); axis equal;
colormap('jet'); caxis([-3e5 2e5]);  colorbar; box on;
title('\sigma_{xy}');
hold off;

%========================== convergence plot ==============================
t = 11;  % First 10 time steps
[~, total_time_steps] = size(error);
if total_time_steps < t
    error("Your 'error' matrix has fewer than 10 time steps.");
end
figure;
hold on;
for j = 1:t
    err_vals = error(:, j+1);
    % Remove zeros (assumes iterations stop after convergence)
    valid_idx = err_vals ~= 0;
    err_vals = err_vals(valid_idx);
    iter_count = numel(err_vals);
    h = plot(j * ones(1, iter_count), err_vals, 'k-o','MarkerFaceColor', 'k', 'MarkerSize', 7);
end
legend(h, 'Iterations', 'Location', 'northeast');
set(gca, 'YScale', 'log');
xlabel('Time steps');
ylabel('e_{disp}');
xlim([0 t+1]); ylim([10e-8 10e0]);
sgtitle('Convergence iterations','Fontsize',10);
grid on; box on;

%{
%=========================== plot for fluid flux ==========================
%__________________________________________________________________________
% Assuming existing variables: elementNodes, nodeCoords, displacement, flux_all, dt
num_elements = size(elementNodes, 1);
element_centers = zeros(num_elements, 2);

% Compute element centers (for flux magnitude interpolation)
for e = 1:num_elements
    node_ids = elementNodes(e, :);  % 8 nodes per Q8 element
    coords = nodeCoords(node_ids, :);  % 8 x 2
    element_centers(e, :) = mean(coords, 1);  % Average of 8 nodes
end

xc = element_centers(:, 1);
yc = element_centers(:, 2);

% Define time steps
dt = 0.01;
time_indices = round([0.1, 0.66, 20] / dt) + 1;  % Adjusted to match figure times

% Create a fine uniform grid for interpolating flux magnitude
x_min = min(nodeCoords(:, 1));
x_max = max(nodeCoords(:, 1));
y_min = min(nodeCoords(:, 2));
y_max = max(nodeCoords(:, 2));
[xq, yq] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));

% Gauss points for 2x2 quadrature in reference element [-1, 1] x [-1, 1]
gaussPoints = [-sqrt(1/3), sqrt(1/3)];  % Gauss points: -0.577, 0.577
numGaussPerElement = 4;  % 2x2 = 4 Gauss points per element

% Loop over time steps
for idx = 1:length(time_indices)
    t_idx = time_indices(idx);
    t = (t_idx - 1) * dt;  % Convert index to actual time

    % Extract flux components for the current time step (already at Gauss points)
    fx_gauss = flux_all(:, :, 1, t_idx);  % x-component at Gauss points (num_elements x 4)
    fy_gauss = flux_all(:, :, 2, t_idx);  % y-component at Gauss points (num_elements x 4)
    
    % Average flux at element centers for contour plot
    fx = squeeze(mean(flux_all(:, :, 1, t_idx), 2));  % x-component per element
    fy = squeeze(mean(flux_all(:, :, 2, t_idx), 2));  % y-component per element
    flux_mag = sqrt(fx.^2 + fy.^2);  % Magnitude per element

    % Interpolate flux magnitude to uniform grid for contour plot
    flux_mag_grid = griddata(xc, yc, flux_mag, xq, yq, 'natural');  % Interpolate magnitude

    % Deformed coordinates for quadrant I
    ux = displacement(1:2:end, t_idx);  % x-displacements at current time
    uy = displacement(2:2:end, t_idx);  % y-displacements at current time
    deformed_coords_I = nodeCoords + [ux, uy];

    % Mirror deformed coordinates for quadrant II
    deformed_coords_II = deformed_coords_I;
    deformed_coords_II(:, 1) = -deformed_coords_I(:, 1);  % quadrant II: (-X, +Y)

    % Compute Gauss point positions in deformed configuration
    gaussPositions_I = zeros(num_elements * numGaussPerElement, 2);  % Deformed positions for quadrant I
    gaussPositions_II = zeros(num_elements * numGaussPerElement, 2); % Deformed positions for quadrant II
    gaussFluxX = zeros(num_elements * numGaussPerElement, 1);  % Flux x-component at Gauss points
    gaussFluxY = zeros(num_elements * numGaussPerElement, 1);  % Flux y-component at Gauss points

    for e = 1:num_elements
        node_ids = elementNodes(e, :);
        coords_deformed_I = deformed_coords_I(node_ids, :);  % Deformed coordinates of element nodes (quadrant I)
        coords_deformed_II = deformed_coords_II(node_ids, :); % Deformed coordinates for quadrant II

        % Loop over Gauss points
        gp_idx = 1;
        for xi = gaussPoints
            for eta = gaussPoints
                % Shape functions for Q8 element at (xi, eta)
                N = zeros(1, 8);
                N(1) = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1);
                N(2) = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1);
                N(3) = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1);
                N(4) = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1);
                N(5) = 0.5 * (1 - xi^2) * (1 - eta);
                N(6) = 0.5 * (1 + xi) * (1 - eta^2);
                N(7) = 0.5 * (1 - xi^2) * (1 + eta);
                N(8) = 0.5 * (1 - xi) * (1 - eta^2);

                % Global position of Gauss point in deformed configuration
                global_idx = (e - 1) * numGaussPerElement + gp_idx;
                gaussPositions_I(global_idx, :) = N * coords_deformed_I;  % Quadrant I
                gaussPositions_II(global_idx, :) = N * coords_deformed_II; % Quadrant II

                % Store flux components at Gauss point
                gaussFluxX(global_idx) = fx_gauss(e, gp_idx);  % x-component
                gaussFluxY(global_idx) = fy_gauss(e, gp_idx);  % y-component

                gp_idx = gp_idx + 1;
            end
        end
    end

    % Create figure
    figure;
    hold on;

    % Plot deformed mesh for quadrant I
    for element_id = 1:min(50, size(elementNodes, 1))
        nodes = elementNodes(element_id, :);
        coords_deformed_I = deformed_coords_I(nodes, :);
        x_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 1);  % Connect nodes in order
        y_def_I = coords_deformed_I([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_I, y_def_I, 'k-', 'LineWidth', 0.1);
    end

    % Plot deformed mesh for quadrant II
    for element_id = 1:min(50, size(elementNodes, 1))
        nodes = elementNodes(element_id, :);
        coords_deformed_II = deformed_coords_II(nodes, :);
        x_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 1);
        y_def_II = coords_deformed_II([1 5 2 6 3 7 4 8 1], 2);
        plot(x_def_II, y_def_II, 'k-', 'LineWidth', 0.1);
    end

    % Overlay flux vectors at Gauss points using quiver
    quiver(gaussPositions_I(:, 1), gaussPositions_I(:, 2), gaussFluxX, gaussFluxY, 0.5, 'r');  % Quadrant I
    quiver(gaussPositions_II(:, 1), gaussPositions_II(:, 2), -gaussFluxX, gaussFluxY, 0.5, 'r');  % Quadrant II (mirror x-component)

    % Labels and title
    title(['Fluid Flux at t = ', num2str(t)]);
    xlabel('x');
    ylabel('y');
    axis equal;  % Ensure aspect ratio matches deformed shape
    hold off;
end
    %}
end
end
%%