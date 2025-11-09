%% 1. 定义FK谱绘制函数（新增纵轴范围参数，支持自定义调整）
function plot_fk_spectrum1(seismic_data, dt, dx, subplot_pos, title_str, freq_lim, x_lim)
    % 计算并在指定子图位置显示地震数据的FK谱（支持纵轴频率范围自定义）
    % 输入:
    %   seismic_data: 地震数据矩阵 (时间采样点 × 道数)
    %   dt: 时间采样间隔 (秒)
    %   dx: 道间距 (米)
    %   subplot_pos: 子图位置（如[2,4,4]表示2行4列第4幅）
    %   title_str: 子图标题
    %   freq_lim: 纵轴频率显示范围，格式为[min_freq, max_freq]（单位：Hz）
    
    [nt, nx] = size(seismic_data);  % nt: 时间采样点数, nx: 道数
    
    % 1. 二维傅里叶变换与幅度转换
    fk_data = fftshift(fft2(seismic_data));
    fk_amp = abs(fk_data);
    fk_db = 20*log10(fk_amp/max(fk_amp(:)));  % 归一化到dB，突出能量分布
    
    % 2. 计算频率轴（Hz）和波数轴（1/m）
    fs = 1/dt;  % 采样频率
    f = (-nt/2:nt/2-1)*(fs/nt);  % 频率轴（含正负频率）
    k = (-nx/2:nx/2-1)*(1/(nx*dx/1000));  % 波数轴（含正负波数）
    
    % 3. 在指定子图位置绘制FK谱
    subplot(subplot_pos(1), subplot_pos(2), subplot_pos(3));
    imagesc(k, f, fk_db);
    axis xy;  % y轴向上为正，符合地震数据显示习惯
%     colormap(gray);
    colorbar;  % 显示颜色条，便于对比能量强弱
    caxis([-60 0]);  % 限制dB范围，过滤微弱噪声，突出有效信号
    
    % 关键：根据输入参数设置纵轴频率范围（支持自定义调整）
    ylim(freq_lim);  
    xlim(x_lim);
    % 标签与格式设置
    xlabel('Wavenumber (1/km)');
    ylabel('Frequency (Hz)');
    title(title_str);
    set(gca,'FontName','Times New Roman','FontSize',10);
end