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

[data_idx, theta] = frame_sync(filtered_rx_signal,conf);

Len = ((conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.repeatTrainingFrequency))+1)*(conf.f_s/conf.f_sym+conf.LengthCP);
RX_Time_Vector = filtered_rx_signal(data_idx:data_idx+Len-1); %the length of our signal is the training and the data 
RX_Time_Matrix = reshape(RX_Time_Vector,conf.f_s/conf.f_sym+conf.LengthCP,(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.repeatTrainingFrequency))+1);


%% channel estimation & phase correction
rx_corrected = channel_estimation(RX_Time_Matrix, conf);

%% demap bits
rxbits = demapper(rx_corrected);


end

 


