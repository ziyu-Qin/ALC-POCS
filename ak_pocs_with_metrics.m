function [result, iter, metrics] = ak_pocs_with_metrics(data, mask, max_iter, ratio, lambda, use_lateral_constraint, diff_type, original_full)
    % Three-convex-set projection interpolation function with iteration metrics recording
    % Records relative error, PSNR, and trace correlation coefficient for each iteration to analyze convergence
    % Inputs:
    %   data - Seismic data with missing traces (zeros in missing positions)
    %   mask - Logical matrix indicating known (true) and missing (false) traces
    %   max_iter - Maximum number of iterations
    %   ratio - Threshold scaling ratio
    %   lambda - Weight for lateral constraint
    %   use_lateral_constraint - Logical flag to enable/disable lateral constraint
    %   diff_type - Type of inter-trace difference matrix (1 or 2)
    %   original_full - Complete original seismic data (for calculating iteration error metrics)
    % Outputs:
    %   result - Interpolation result
    %   iter - Actual number of iterations performed
    %   metrics - Structure containing iteration metrics:
    %       rel_error: Relative error per iteration (Frobenius norm ratio of error to original data)
    %       psnr: Peak signal-to-noise ratio per iteration (dB)
    %       trace_cc: Mean correlation coefficient between interpolated and original data (per trace)
    
    % Copyright (c) [2025] [Ziyu Qin]
    % All rights reserved.
    % For academic/research use only.  
    % -------------------------- 1. Parameter initialization --------------------------
    if nargin < 6
        use_lateral_constraint = true;  % Default: enable lateral constraint
    end
    if nargin < 7
        diff_type = 1; % Default: use difference matrix type 1
    end
    if nargin < 8
        error('Input original_full is required for calculating iteration metrics!');
    end
    
    [cols, rows] = size(data);  % cols: time samples; rows: number of traces
    result = data;
    result(~mask) = 0;  % Initialize missing values to zero
    iter = 0;
    
    % Build inter-trace difference matrix D
    if use_lateral_constraint
        if diff_type == 1
            D = zeros(rows, rows);
            for i = 2 : rows - 1
                D(i, i - 1) = -1;
                D(i, i) = 2;
                D(i, i + 1) = -1;
            end
            D(1, 1) = 1;
            D(1, 2) = -1;
            D(end, end - 1) = -1;
            D(end, end) = 1;
        elseif diff_type == 2
            D = zeros(rows, rows);
            for i = 1 : rows - 1
                D(i, i) = 1;
                D(i, i + 1) = -1;
            end
            D(end, end) = 1;
        else
            error('Invalid diff_type. Choose 1 or 2.');
        end
    end
    
    % Calculate exponentially decaying threshold sequence for AK-POCS
    initial_A = project_A(result, data, mask);
    fft_amp = abs(fftshift(fft2(initial_A)));
    p_max = max(fft_amp(:)) * ratio;
    p_min = p_max * 1e-4;
    
    if max_iter == 1
        thresholds = p_max;
    else
        b = -log(p_max / p_min) / (max_iter - 1);
        thresholds = p_max * exp(b * (0:max_iter-1)');
    end
    
    % -------------------------- 2. Initialize iteration metric arrays --------------------------
    metrics.rel_error = zeros(1, max_iter);  % Relative error (preallocated for max iterations)
    metrics.psnr = zeros(1, max_iter);       % PSNR (dB)
    metrics.trace_cc = zeros(1, max_iter);   % Mean trace correlation coefficient
    data_range = max(original_full(:)) - min(original_full(:));  % Data amplitude range for PSNR
    if data_range == 0
        data_range = 1;  % Avoid division by zero in PSNR calculation
    end
    
    % -------------------------- 3. Main iteration loop --------------------------
    while iter < max_iter
        prev = result;
        iter = iter + 1;
        t = thresholds(iter);  % Threshold for current iteration
        
        % Frequency domain constraint
        result = project_B(result, t);
        
        % Lateral constraint (if enabled)
        if use_lateral_constraint
            diff = result * D;  % Calculate inter-trace differences
            grad_diff = sign(diff) * D';  % Subgradient of L1 norm
            result = result - 0.001 * lambda * grad_diff;  % Gradient update
        end
 
        % Convergence check
        if norm(result - prev, 'fro') < 1e-6
            break;
        end
        
        % Known data constraint (preserve known values)
        result = project_A(result, data, mask);
        
        % -------------------------- Calculate current iteration metrics --------------------------
        % Relative error
        error_mat = result - original_full;
        metrics.rel_error(iter) = norm(error_mat, 'fro') / norm(original_full, 'fro');
        
        % PSNR
        mse = mean(error_mat(:).^2);  % Mean squared error
        metrics.psnr(iter) = 10 * log10((data_range^2) / mse);
        
        % Mean trace correlation coefficient (exclude NaN values)
        cc_values = arrayfun(@(j) corr(result(:,j), original_full(:,j)), 1:rows);
        metrics.trace_cc(iter) = mean(cc_values(~isnan(cc_values)));
        
        % Print iteration information
        if mod(iter, 5) == 0 || iter == max_iter
            if use_lateral_constraint
                fprintf('Iteration %d/%d (with lateral constraint) | Rel Error: %.6f | PSNR: %.2f dB | Trace CC: %.4f\n', ...
                    iter, max_iter, metrics.rel_error(iter), metrics.psnr(iter), metrics.trace_cc(iter));
            else
                fprintf('Iteration %d/%d (no lateral constraint) | Rel Error: %.6f | PSNR: %.2f dB | Trace CC: %.4f\n', ...
                    iter, max_iter, metrics.rel_error(iter), metrics.psnr(iter), metrics.trace_cc(iter));
            end
        end
    end
    
    % -------------------------- Truncate metrics to actual iterations --------------------------
    metrics.rel_error = metrics.rel_error(1:iter);
    metrics.psnr = metrics.psnr(1:iter);
    metrics.trace_cc = metrics.trace_cc(1:iter);
    
    % -------------------------- Internal subfunctions --------------------------
    % Data domain projection: preserve known data values
    function proj_A = project_A(current, original, mask)
        proj_A = current;
        proj_A(mask) = original(mask);
    end
    
    % Frequency domain projection: apply thresholding in frequency domain
    function proj_B = project_B(current, threshold)
        fft_data = fftshift(fft2(current));
        fft_data(abs(fft_data) < threshold) = 0;
        proj_B = real(ifft2(ifftshift(fft_data)));
    end
end