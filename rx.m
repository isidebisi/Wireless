function [rxbits rx_corrected conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

if strcmp(conf.plotfigure,'true')
    fftSignal = abs(fftshift(fft(rxsignal)));
    N = length(fftSignal); 
    f = (-N/2:N/2-1) * (conf.f_s/N);
    
    figure(7);
    plot(f,fftSignal);
    title('Spectrum of transmitted signal after modulation');
    xlabel('frequency (Hz)');
    ylabel('Amplitude');
    %xlim([-10000 10000])
end   

%% demodulate
demodulated_signal = demodulate(rxsignal, conf);

%% low pass filter
filtered_rx_signal = ofdmlowpass(demodulated_signal,conf,conf.enlarged_bandwidth);

%% Frame synchronization
% retrieve start of data index
[data_index, theta] = frame_sync(filtered_rx_signal, conf);

% Remove preamble
signal_length = ((conf.OFDM_symbols + floor(conf.OFDM_symbols / conf.f_train)) + 1) * (conf.f_s / conf.spacing + conf.cp_len);

% extract signal
received_signal = filtered_rx_signal(data_index:data_index + signal_length - 1); 

%% remove cyclic prefix
% Reshape the received signal into a matrix
time_matrix = reshape(received_signal, conf.f_s / conf.spacing + conf.cp_len, (conf.OFDM_symbols + floor(conf.OFDM_symbols / conf.f_train)) + 1);

% remove cp
rx_no_cp = time_matrix(conf.cp_len + 1:end, :);

%% channel estimation, phase correction & frequency domain conversion
rx_corrected = channel_estimation(rx_no_cp, conf);

%% demapper QPSK
rxbits = demapper(rx_corrected);


end

 


