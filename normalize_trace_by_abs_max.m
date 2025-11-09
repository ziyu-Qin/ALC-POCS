function normalized_data = normalize_trace_by_abs_max(data)
    % 地震道归一化：每道除以该道绝对值的最大值，使每道最大绝对值为1
    % 输入：
    %   data - 原始地震数据矩阵，维度为 [采样点数量, 道数]（行=采样点，列=地震道）
    % 输出：
    %   normalized_data - 归一化后的数据矩阵，维度与输入完全一致
    
    % 1. 获取数据维度与初始化输出矩阵
    [samples, traces] = size(data);
    normalized_data = zeros(size(data));
    
    % 2. 逐道进行归一化处理
    for trace_idx = 1:traces
        % 提取当前地震道数据
        current_trace = data(:, trace_idx);
        
        % 计算当前道绝对值的最大值（避免NaN影响）
        abs_max = max(abs(current_trace), [], 'omitnan');
        
        % 3. 避免除零错误（若道数据全为0，保持全0）
        if abs_max < 1e-10
            normalized_data(:, trace_idx) = current_trace;
        else
            % 核心归一化：当前道除以自身绝对值最大值
            normalized_data(:, trace_idx) = current_trace / abs_max;
        end
    end
end