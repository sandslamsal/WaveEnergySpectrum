function [Hm0, f, S, Tp, varargout] = spectrumWaves(varargin)

% [spectra, frequencies, Hm0, Tp] = spectrumWaves(data1, data2, ..., data6, delt);
% Author: Sandesh Lamsal
% Description: Compute energy spectra and plot time series and spectra of wave height data up to 6 different Wave Data.


    % Check for potential errors in the inputs
    num_datasets = length(varargin) - 1; % Exclude delt from the count
    if num_datasets < 1 || num_datasets > 6
        error('The function accepts between 1 and 6 datasets.');
    end
    
    delt = varargin{end}; % Last argument is delt
    
    for i = 1:num_datasets
        if isempty(varargin{i})
            error('Input datasets must not be empty.');
        end
    end
    
    if delt <= 0
        error('Input delt must be positive.');
    end

    N = length(varargin{1});
    
    % Compute Fourier coefficients and spectra for each dataset
    S = cell(1, num_datasets);
    f = cell(1, num_datasets);
    Hm0 = zeros(1, num_datasets);
    Tp = zeros(1, num_datasets);
    
    for i = 1:num_datasets
        z = varargin{i};
        zmean = mean(z, 'omitnan');
        zsub = z - zmean;
        
        Ce = fft(zsub, N);
        S{i} = Ce.*conj(Ce) / N;
        S{i}(round(N/2) + 2:end) = [];
        delf = 1 / (delt * N);
        S{i}(2:round(N/2)) = 2 * S{i}(2:round(N/2));
        f{i} = delf * (0:length(S{i}) - 1);
        
        % Calculate variance and Hm0
        var = sum(S{i}) * delf;
        Hm0(i) = 4.004 * sqrt(var);
        [~, ind] = max(S{i});
        fmax = f{i}(ind);
        Tp(i) = 1 / fmax;
    end
    
    % Output results
    varargout = [S, f];
    
    % Plot time series
    figure;
    for i = 1:num_datasets
        subplot(num_datasets, 1, i);
        plot((0:delt:(N - 1) * delt), varargin{i}, 'LineWidth', 1.5);
        title(['Time Series - Dataset ', num2str(i)]);
        xlabel('Time (s)');
        ylabel('{\eta} (m)');
        legend(['Wave Probe', num2str(i)]);
        grid on;
    end
    
    % Plot energy spectrum for each dataset
    figure;
    for i = 1:num_datasets
        indices = f{i} <= 2;
        subplot(num_datasets, 1, i);
        plot(f{i}(indices), S{i}(indices), 'LineWidth', 1.5);
        title(['Energy Spectrum - Dataset ', num2str(i)]);
        xlabel('Frequency (Hz)');
        ylabel('S_{\eta} (m^2/Hz)');
        legend(['Wave Probe', num2str(i)]);
        grid on;
    end
    %{
    % Plot combined energy spectrum
    figure;
    hold on;
    for i = 1:num_datasets
        indices = f{i} <= 2;
        plot(f{i}(indices), S{i}(indices), 'LineWidth', 1.5);
    end
    hold off;
    title('Combined Energy Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('S_{\eta} (m^2/Hz)');
    legend('WP 1', 'WP 2', 'WP 3', 'WP 4', 'WP 5', 'WP 6');
    grid on;
%}
end
