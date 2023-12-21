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

% dummy 

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

rxbits = zeros(conf.nbits,1);

%% demodulate
demodulated_signal = demodulate(rxsignal, conf);

%% low pass filter
filtered_rx_signal = ofdmlowpass(demodulated_signal,conf,conf.enlarged_bandwidth);

%% Frame synchronization

% start of frame detection 
[data_idx, theta] = frame_sync(filtered_rx_signal,conf);

% remove preamble
Len = ((conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))+1)*(conf.f_s/conf.f_sym+conf.cp_len);
RX_Time_Vector = filtered_rx_signal(data_idx:data_idx+Len-1); %the length of our signal is the training and the data 
RX_Time_Matrix = reshape(RX_Time_Vector,conf.f_s/conf.f_sym+conf.cp_len,(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))+1);

%% remove cyclic prefix

rx_no_cp = RX_Time_Matrix(conf.cp_len + 1:end, :);

%% channel estimation & phase correction
rx_corrected = channel_estimation(rx_no_cp, conf);

%% demapper QPSK
rxbits = demapper(rx_corrected);


end

 


