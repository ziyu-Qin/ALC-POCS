% Interpolation Comparison Script (With Iteration Metrics)
% Copyright (c) [2025] [Ziyu Qin]
% All rights reserved.
% For academic/research use only.  

clear; clc; close all;

%% 1. Load and preprocess synthetic seismic data
load('theoretical_pre_data.mat');  
synthetic = theoretical_pre_data;
[rows, cols] = size(synthetic);  % rows=time samples, cols=number of traces
[cmin, cmax] = deal(min(synthetic(:)), max(synthetic(:)));  % Data amplitude range

%% 2. Key parameter settings (1ms sampling)
sample_rate_ms = 1;  % Time sampling rate (1ms/point)
dt = sample_rate_ms / 1000;  % Time interval (s): 1ms = 0.001s
dx = 25;  % Trace spacing (m, adjust for real data)
t_ms = (1:rows) * sample_rate_ms;  % Time axis (ms)
trace_num = 1:cols;  % Trace number axis

% FK spectrum vertical frequency range (adjustable)
freq_range_akpocs = [-50, 50];  
freq_range_this = [-50, 50];  
x_range=[-20,20];
% Print basic data info
fprintf('Data loaded successfully!\n');
fprintf('  - Time sample rate: %.0f ms/point (dt = %.3f s)\n', sample_rate_ms, dt);
fprintf('  - Trace spacing: %.0f m\n', dx);
fprintf('  - Data size: Time samples = %d, Traces = %d\n', rows, cols);
fprintf('  - FK freq range: %.0f ~ %.0f Hz\n', freq_range_akpocs(1), freq_range_akpocs(2));

%% 3. Simulate data missing (30% missing)
missing_ratio = 0.3;  % 30% traces missing
rng(6);  % Fix random seed for reproducibility
missing_traces = randperm(cols, round(missing_ratio * cols));  % Randomly select missing traces

% Generate missing data and mask
mask = true(rows, cols);
mask(:, missing_traces) = false;
data_missing = synthetic;
data_missing(~mask) = 0;

%% 4. Interpolation algorithm parameters
max_iter = 100;       % Maximum iterations
ratio = 1;            % Threshold ratio
tv_lambda = 5;        % Inter-trace constraint weight

% AK-POCS (without inter-trace constraint): Call external function with metrics
[result_akpocs, iter_akpocs, metrics_akpocs] = ak_pocs_with_metrics(...
    data_missing, mask, max_iter, ratio, tv_lambda, false, 1, synthetic);
 
[result_this, iter_this, metrics_this] = alc_pocs_with_metrics(data_missing, mask, max_iter,...
    ratio, 50, 50, 0.70, synthetic);

%% 6. Result visualization (8 interpolation plots + 3 iteration metric plots)
% -------------------------- Figure 1: 8 interpolation result plots --------------------------
figure('Position', [100, 100, 1800, 600]);
set(gcf, 'Color', 'white');

% Top row subplots (1-4)
% Subplot 1: Original data
% subplot(2,4,1);
% imagesc(trace_num, t_ms, synthetic);
% caxis([cmin, cmax]);
% colormap(gray);
% title('Original Data');
% xlabel('Trace Number');
% ylabel('Time (ms)');
% set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 2: AK-POCS interpolation result
subplot(2,4,1);
imagesc(trace_num, t_ms, result_akpocs);
caxis([cmin, cmax]);
colormap(gray);
title('AK-POCS Interpolation');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 3: AK-POCS interpolation error
subplot(2,4,2);
err_akpocs = result_akpocs - synthetic;
imagesc(trace_num, t_ms, err_akpocs);
colormap(gray);
caxis([-max(abs(err_akpocs(:))), max(abs(err_akpocs(:)))]);  % Symmetric error display
title('AK-POCS Error');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 4: AK-POCS amplified error
subplot(2,4,3);
imagesc(trace_num, t_ms, err_akpocs);
colormap(gray);
caxis([-max(abs(err_akpocs(:))/3), max(abs(err_akpocs(:)))/3]);  % Amplified symmetric error
title('AK-POCS Amplified Error');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 5: AK-POCS FK spectrum (call existing plot_fk_spectrum1)
plot_fk_spectrum1(result_akpocs, dt, dx, [2,4,4], 'AK-POCS FK Spectrum', freq_range_akpocs, x_range);

% Bottom row subplots (5-8)
% Subplot 6: Missing data
% subplot(2,4,5);
% imagesc(trace_num, t_ms, data_missing);
% colormap(gray);
% caxis([cmin, cmax]);
% title(sprintf('Missing Data (%.0f%% Missing)', missing_ratio*100));
% xlabel('Trace Number');
% ylabel('Time (ms)');
% set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 7: Proposed method interpolation result
subplot(2,4,5);
imagesc(trace_num, t_ms, result_this);
caxis([cmin, cmax]);
colormap(gray);
title('Proposed Method Interpolation');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 8: Proposed method interpolation error
subplot(2,4,6);
err_this = result_this - synthetic;
imagesc(trace_num, t_ms, err_this);
colormap(gray);
caxis([-max(abs(err_this(:))), max(abs(err_this(:)))]);  % Symmetric error display
title('Proposed Method Error');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 9: Proposed method amplified error
subplot(2,4,7);
imagesc(trace_num, t_ms, err_this);
colormap(gray);
caxis([-max(abs(err_this(:))/3), max(abs(err_this(:))/3)]);  % Amplified symmetric error
title('Proposed Method Amplified Error');
xlabel('Trace Number');
ylabel('Time (ms)');
set(gca,'FontName','Times New Roman','FontSize',10); 

% Subplot 10: Proposed method FK spectrum (call existing plot_fk_spectrum1)
plot_fk_spectrum1(result_this, dt, dx, [2,4,8], 'Proposed Method FK Spectrum', freq_range_this, x_range);

% -------------------------- Figure 2: Iteration metric curves --------------------------
figure('Position', [200, 200, 1600, 380]);
set(gcf, 'Color', 'white');

% Unified iteration vector (align x-axis for two curves)
iter_vec_akpocs = 1:iter_akpocs;
iter_vec_this = 1:iter_this;

% Subplot 1: Relative error vs iteration
subplot(1,3,1);
plot(iter_vec_akpocs, metrics_akpocs.rel_error, 'b-', 'LineWidth', 1.5, 'DisplayName', 'AK-POCS');
hold on;
plot(iter_vec_this, metrics_this.rel_error, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Proposed Method');
xlabel('Iteration Number');
ylabel('Relative Error');
title('Relative Error vs Iteration');
legend('Location', 'Best');
grid on;
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 2: PSNR vs iteration
subplot(1,3,2);
plot(iter_vec_akpocs, metrics_akpocs.psnr, 'b-', 'LineWidth', 1.5, 'DisplayName', 'AK-POCS');
hold on;
plot(iter_vec_this, metrics_this.psnr, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Proposed Method');
xlabel('Iteration Number');
ylabel('PSNR (dB)');
title('PSNR vs Iteration');
legend('Location', 'Best');
grid on;
set(gca,'FontName','Times New Roman','FontSize',10);

% Subplot 3: Trace cross-correlation vs iteration
subplot(1,3,3);
plot(iter_vec_akpocs, metrics_akpocs.trace_cc, 'b-', 'LineWidth', 1.5, 'DisplayName', 'AK-POCS');
hold on;
plot(iter_vec_this, metrics_this.trace_cc, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Proposed Method');
xlabel('Iteration Number');
ylabel('Trace Cross-Correlation');
title('Trace CC vs Iteration');
legend('Location', 'Best');
grid on;
set(gca,'FontName','Times New Roman','FontSize',10);

%% 7. Quantitative evaluation calculation and print (final results)
% Relative error
error_this = norm(result_this - synthetic, 'fro') / norm(synthetic, 'fro');
error_akpocs = norm(result_akpocs - synthetic, 'fro') / norm(synthetic, 'fro');

% PSNR (avoid division by zero)
data_range = max(synthetic(:)) - min(synthetic(:));
if data_range == 0
    data_range = 1;
end
psnr_this = 10 * log10((rows*cols*data_range^2) / norm(result_this - synthetic, 'fro')^2);
psnr_akpocs = 10 * log10((rows*cols*data_range^2) / norm(result_akpocs - synthetic, 'fro')^2);

% Trace cross-correlation (exclude NaN)
cc_values_this = arrayfun(@(j) corr(result_this(:,j), result_this(:,j+1)), 1:cols-1);
cc_values_akpocs = arrayfun(@(j) corr(result_akpocs(:,j), result_akpocs(:,j+1)), 1:cols-1);
cc_this = mean(cc_values_this(~isnan(cc_values_this)));
cc_akpocs = mean(cc_values_akpocs(~isnan(cc_values_akpocs)));

% Print results
fprintf('\n=========================================\n');
fprintf('Interpolation Evaluation Results \n');
fprintf('=========================================\n');
fprintf('| Method         | Rel Error | PSNR (dB) | Trace CC |\n');
fprintf('|----------------|-----------|-----------|----------|\n');
fprintf('| Proposed       | %.6f     | %.2f      | %.3f     |\n', error_this, psnr_this, cc_this);
fprintf('| AK-POCS        | %.6f     | %.2f      | %.3f     |\n', error_akpocs, psnr_akpocs, cc_akpocs);

fprintf('=========================================\n');
