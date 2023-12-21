function [rxbits conf] = rx(rxsignal,conf,k)

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

%% Demodulate
demodulated_signal = demodulate(rxsignal, conf);

%% Low pass filter
filtered_rx_signal = ofdmlowpass(demodulated_signal,conf,conf.enlarged_bandwidth);

%% Frame synchronization
% retrieve start of data index
[data_index, theta] = frame_sync(filtered_rx_signal, conf);

% Remove preamble
signal_length = ((conf.ofdm_symbols + floor(conf.ofdm_symbols / conf.f_train)) + 1) * (conf.f_s / conf.f_sym + conf.cp_len);

% extract signal
received_signal = filtered_rx_signal(data_index:data_index + signal_length - 1); 

%% Cyclic prefix removal
% Reshape the received signal into a matrix
time_matrix = reshape(received_signal, conf.f_s / conf.f_sym + conf.cp_len, (conf.ofdm_symbols + floor(conf.ofdm_symbols / conf.f_train)) + 1);

% remove cp
rx_no_cp = time_matrix(conf.cp_len + 1:end, :);

%% Channel estimation, phase correction & frequency domain conversion
rx_corrected = channel_estimation(rx_no_cp, conf);

%% Demapper QPSK
rxbits = demapper(rx_corrected);

end

 


