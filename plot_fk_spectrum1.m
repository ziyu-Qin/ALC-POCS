%% 1. Define FK spectrum plotting function (add vertical axis range parameter for custom adjustment)
function plot_fk_spectrum1(seismic_data, dt, dx, subplot_pos, title_str, freq_lim, x_lim)
    % Calculate and display FK spectrum of seismic data at specified subplot position (supports custom vertical frequency range)
    % Inputs:
    %   seismic_data: Seismic data matrix (time samples × number of traces)
    %   dt: Time sampling interval (seconds)
    %   dx: Trace spacing (meters)
    %   subplot_pos: Subplot position (e.g., [2,4,4] = 2 rows, 4 columns, 4th subplot)
    %   title_str: Subplot title
    %   freq_lim: Vertical frequency display range, format [min_freq, max_freq] (unit: Hz)
    %   x_lim: Horizontal wavenumber display range, format [min_k, max_k] (unit: 1/km)
    
    [nt, nx] = size(seismic_data);  % nt: number of time samples, nx: number of traces
    
    % 1. 2D Fourier transform and amplitude conversion
    fk_data = fftshift(fft2(seismic_data));
    fk_amp = abs(fk_data);
    fk_db = 20*log10(fk_amp/max(fk_amp(:)));  % Normalize to dB to highlight energy distribution
    
    % 2. Calculate frequency axis (Hz) and wavenumber axis (1/km)
    fs = 1/dt;  % Sampling frequency
    f = (-nt/2:nt/2-1)*(fs/nt);  % Frequency axis (includes positive/negative frequencies)
    k = (-nx/2:nx/2-1)*(1/(nx*dx/1000));  % Wavenumber axis (includes positive/negative wavenumbers)
    
    % 3. Plot FK spectrum at specified subplot position
    subplot(subplot_pos(1), subplot_pos(2), subplot_pos(3));
    imagesc(k, f, fk_db);
    axis xy;  % Invert y-axis to align with seismic data display convention
%     colormap(gray);
    colorbar;  % Show colorbar for energy intensity comparison
    caxis([-60 0]);  % Limit dB range to filter weak noise and emphasize valid signals
    
    % Key: Set vertical frequency range via input parameter (supports customization)
    ylim(freq_lim);  
    xlim(x_lim);
    % Label and format settings
    xlabel('Wavenumber (1/km)');
    ylabel('Frequency (Hz)');
    title(title_str);
    set(gca,'FontName','Times New Roman','FontSize',10);
end
