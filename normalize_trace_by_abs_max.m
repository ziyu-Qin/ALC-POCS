function normalized_data = normalize_trace_by_abs_max(data)
    % Normalize seismic traces by dividing each trace by its absolute maximum value, making the max absolute value of each trace 1
    % Input:
    %   data - Original seismic data matrix with dimensions [number of samples, number of traces] (rows = samples, columns = traces)
    % Output:
    %   normalized_data - Normalized data matrix with the same dimensions as input
    
    % 1. Get data dimensions and initialize output matrix
    [samples, traces] = size(data);
    normalized_data = zeros(size(data));
    
    % 2. Normalize each trace individually
    for trace_idx = 1:traces
        % Extract current trace data
        current_trace = data(:, trace_idx);
        
        % Calculate the maximum absolute value of the current trace (ignore NaN)
        abs_max = max(abs(current_trace), [], 'omitnan');
        
        % 3. Avoid division by zero (keep all zeros if the trace is all zeros)
        if abs_max < 1e-10
            normalized_data(:, trace_idx) = current_trace;
        else
            % Core normalization: divide current trace by its absolute maximum
            normalized_data(:, trace_idx) = current_trace / abs_max;
        end
    end
end
