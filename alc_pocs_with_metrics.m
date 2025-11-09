function [result, iter, metrics] = alc_pocs_with_metrics(data, mask, max_iter, ratio, lambda, win_size, sim_threshold, original_full)
    % ALC-POCS algorithm with iteration metrics recording
    % Records relative error, PSNR, and trace correlation coefficient for each iteration to analyze convergence
    % Inputs:
    %   data - Seismic data with missing traces (zeros in missing positions)
    %   mask - Logical matrix: true for known traces, false for missing traces
    %   max_iter - Maximum number of iterations
    %   ratio - Threshold scaling ratio
    %   lambda - Weight for lateral constraint
    %   win_size - Size of non-overlapping time windows
    %   sim_threshold - Threshold for valid similarity values
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
    % -------------------------- 1. Parameter check and default values --------------------------
    if nargin < 6, win_size = 11; end
    if nargin < 7, sim_threshold = 0.6; end
    if nargin < 8
        error('Input original_full is required for calculating iteration metrics!');
    end
    
    [samples, traces] = size(data);  % samples (rows) ¡Á traces (columns)
    result = data;
    result(~mask) = 0;  % Initialize missing values to zero
    iter = 0;
    known_traces = find(mask(1, :));  % Indices of known traces
    if length(known_traces) < 2, error('At least 2 known traces required to find nearest neighbors'); end
    
    % -------------------------- 2. Initialize iteration metric arrays --------------------------
    metrics.rel_error = zeros(1, max_iter);  % Relative error (Frobenius norm ratio)
    metrics.psnr = zeros(1, max_iter);       % Peak signal-to-noise ratio (dB)
    metrics.trace_cc = zeros(1, max_iter);   % Mean trace correlation (interpolation vs original)
    % Data range for PSNR calculation (avoid division by zero)
    data_range = max(original_full(:)) - min(original_full(:));
    if data_range == 0, data_range = 1; end
    
    % -------------------------- 3. Core logic: Inter-trace matrix, time windows, and similarity matrix --------------------------
    % Build inter-trace difference matrix D (traces¡Átraces)
    D = zeros(traces, traces);
    for i = 2 : traces - 1, D(i, [i-1,i,i+1]) = [-1,2,-1]; end
    D(1,1:2) = [1,-1]; D(end,end-1:end) = [-1,1];
    
    % Precompute non-overlapping time window indices and similarity matrix
    [win_start_list, win_end_list, win_num] = gen_nonoverlap_windows(samples, win_size);
    normalized_data = normalize_trace_by_abs_max(data);
    sim_matrix = precompute_similarity_nonoverlap(normalized_data, samples, win_num, win_start_list, win_end_list, traces, known_traces);
    sim_matrix(isnan(sim_matrix)) = 1;  % Handle NaN values
    sim_matrix(sim_matrix == 0) = 1;    % Avoid zero value influence
    
    % Gaussian smoothing for similarity matrix (optional)
    Q = fspecial('gaussian', [20,20], 1e3); 
    sim_matrix = imfilter(sim_matrix, Q, 'replicate'); 
    figure;imagesc(sim_matrix);
    xlabel('Trace Number');
    ylabel('Time (ms)');
    % Generate frequency domain threshold sequence
    initial_A = project_A(result, data, mask);
    fft_amp = abs(fftshift(fft2(initial_A)));
    p_max = max(fft_amp(:)) * ratio;
    p_min = p_max * 1e-4;
    thresholds = p_max * exp(-log(p_max/p_min)/(max_iter-1) * (0:max_iter-1)');

    % -------------------------- 4. Main iteration loop --------------------------
    while iter < max_iter
        prev = result;
        iter = iter + 1;
        t = thresholds(iter);

        % Apply frequency domain constraint
        result = project_B(result, t);

        % Apply dynamic lateral constraint
        if iter > 0
            constraint_matrix = build_constraint_nonoverlap(...
                result, sim_matrix, sim_threshold, win_start_list, win_end_list, win_num, traces);
            result = result - 0.001 * lambda * (result - constraint_matrix) * D';
        end

        % Convergence check and data fidelity constraint
        if norm(result - prev, 'fro') < 1e-6, break; end
        result = project_A(result, data, mask);

        % -------------------------- Calculate current iteration metrics --------------------------
        % Relative error: norm of error divided by norm of original data
        error_mat = result - original_full;
        metrics.rel_error(iter) = norm(error_mat, 'fro') / norm(original_full, 'fro');
        
        % PSNR: calculated based on data amplitude range
        mse = mean(error_mat(:).^2);  % Mean squared error
        metrics.psnr(iter) = 10 * log10((data_range^2) / mse);
        
        % Trace correlation: mean of per-trace correlations (exclude NaN values)
        cc_values = arrayfun(@(j) corr(result(:,j), original_full(:,j)), 1:traces);
        metrics.trace_cc(iter) = mean(cc_values(~isnan(cc_values)));

        % Print iteration information
        if mod(iter, 5) == 0 || iter == max_iter
            fprintf('Iteration %d/%d | Rel Error: %.6f | PSNR: %.2f dB | Trace CC: %.4f\n', ...
                iter, max_iter, metrics.rel_error(iter), metrics.psnr(iter), metrics.trace_cc(iter));
        end
    end

    % -------------------------- Truncate metrics to actual number of iterations --------------------------
    metrics.rel_error = metrics.rel_error(1:iter);
    metrics.psnr = metrics.psnr(1:iter);
    metrics.trace_cc = metrics.trace_cc(1:iter);

end

% -------------------------- Helper functions --------------------------
function [win_start_list, win_end_list, win_num] = gen_nonoverlap_windows(total_samples, win_size)
    % Generate indices for non-overlapping time windows
    win_start_list = 1 : win_size : total_samples;
    win_num = length(win_start_list);
    win_end_list = win_start_list + win_size - 1;
    if win_end_list(end) > total_samples
        win_end_list(end) = total_samples;
    end
end

function sim_matrix = precompute_similarity_nonoverlap(data, samples, win_num, win_start, win_end, traces, known_traces)
    % Precompute similarity matrix using non-overlapping time windows
    sim_matrix = zeros(samples, traces);  % Size: samples¡Átraces (matches input data)

    for trace = 1:traces  % Iterate over each trace (column)
        % Find two nearest known traces for current trace
        dist = abs(known_traces - trace);
        [~, idx] = sort(dist);
        r1 = known_traces(idx(1));
        r2 = known_traces(idx(2));

        % Compute similarity for each window and assign to all samples in the window
        for win_idx = 1:win_num
            s = win_start(win_idx);  % Window start sample (row)
            e = win_end(win_idx);    % Window end sample (row)
            data_r1 = data(s:e, r1);  % Data of known trace r1 in current window
            data_r2 = data(s:e, r2);  % Data of known trace r2 in current window

            % Calculate similarity (single value) for current window
            if length(data_r1) >= 2 && length(data_r2) >= 2
                mean1 = mean(data_r1);
                mean2 = mean(data_r2);
                covar = sum((data_r1-mean1).*(data_r2-mean2))/(length(data_r1)-1);
                std1 = std(data_r1);
                std2 = std(data_r2);
                if std1 > 1e-10 && std2 > 1e-10
                    corr_val = covar/(std1*std2);  % Correlation coefficient (range: [-1,1])
                    corr_val = abs(corr_val);  % Convert to similarity (range: [0,1])
                else
                    corr_val = 0;
                end
            else
                corr_val = 0;
            end

            % Assign window similarity to all samples in the window
            sim_matrix(s:e, trace) = corr_val;
        end
    end
end

function constraint_matrix = build_constraint_nonoverlap(current, sim_matrix, sim_threshold, win_start, win_end, win_num, traces)
    % Build constraint matrix using non-overlapping time windows
    constraint_matrix = current;  % Initialize with current iteration result
    valid_mask = (sim_matrix >= sim_threshold) & ~isnan(sim_matrix);  % Valid sample mask (samples¡Átraces)

    for trace = 1:traces  % Iterate over each trace (column)
        % Find valid samples for current trace
        valid_samples = find(valid_mask(:, trace));
        if isempty(valid_samples), continue; end

        % Determine adjacent traces for constraint
        left = trace - 1;
        right = trace + 1;
        if left < 1
            target = right;
        elseif right > traces
            target = left;
        else
            target = [left, right];
        end

        % Apply constraint to windows with valid samples using adjacent trace data
        for win_idx = 1:win_num
            s = win_start(win_idx);
            e = win_end(win_idx);
            % Check if current window contains valid samples
            if any(valid_samples >= s & valid_samples <= e)
                % Fill window with adjacent trace data
                if length(target) == 2
                    constraint_matrix(s:e, trace) = mean(current(s:e, target), 2);
                else
                    constraint_matrix(s:e, trace) = current(s:e, target);
                end
            end
        end
    end
end

function proj_A = project_A(current, original, mask)
    % Data fidelity projection: retain known data values
    proj_A = current;
    proj_A(mask) = original(mask);
end

function proj_B = project_B(current, threshold)
    % Frequency domain thresholding projection
    fft_data = fftshift(fft2(current));
    fft_data(abs(fft_data) < threshold) = 0;
    proj_B = real(ifft2(ifftshift(fft_data)));
end