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


%% Channel estimation and initial phase estimation using the training data
TrainingData = preamble_generate(conf.N);
TrainingData = -2*(TrainingData) + 1;


payload_data = zeros(conf.N*conf.OFDM_symbols,1);
k =1; % k is the OFDM symbol index which represents a column in the matrix
i =1; % i is the QPSK symbol index which represents the index in a vector 

if strcmp(conf.estimationtype,'block')
    while i < length(payload_data)
        if mod(k,conf.repeatTrainingFrequency+1) == 1
            Y = osfft(RX_Time_Matrix(conf.LengthCP+1:end,k),conf.os_factor);
            H = Y./TrainingData;
            AngleVector = angle(H);
            k= k+1;
        else
        
        payload_data(i:i+conf.N-1) = osfft(RX_Time_Matrix(conf.LengthCP+1:end,k),conf.os_factor);
        payload_data(i:i+length(abs(H))-1) = payload_data(i:i+length(abs(H))-1).*exp(-1j*AngleVector)./abs(H);
        k = k+1;
        i = i + conf.N;
        end
    end

elseif strcmp(conf.estimationtype,'none')
     while i < length(payload_data)
        if mod(k,conf.repeatTrainingFrequency+1) == 1
            k= k+1;
        else
        payload_data(i:i+conf.N-1) = osfft(RX_Time_Matrix(conf.LengthCP+1:end,k),conf.os_factor);
        k = k+1;
        i = i + conf.N;
        end
    end



elseif strcmp(conf.estimationtype,'viterbi')

        while i < length(payload_data)
        if mod(k,conf.repeatTrainingFrequency+1) == 1
            Y = osfft(RX_Time_Matrix(conf.LengthCP+1:end,k),conf.os_factor);
            H = Y./TrainingData;
            AngleVector = angle(H);
            PreviousAngle = angle(H);
            k= k+1;
        else
        
        payload_data(i:i+conf.N-1) = osfft(RX_Time_Matrix(conf.LengthCP+1:end,k),conf.os_factor);
        Segment = payload_data(i:i+conf.N-1);
        for j = 1:conf.N
           A = pi/2*(-1:4);
           deltaTheta = 1/4*angle(Segment(j)^4) + A;
           [~, ind] = min(abs(deltaTheta - PreviousAngle(j)));
           AngleVector(j) = deltaTheta(ind);
           AngleVector(j) = mod(0.01*AngleVector(j) + 0.99*PreviousAngle(j), 2*pi);
           PreviousAngle(j) = AngleVector(j);
        end
        
    
        payload_data(i:i+length(abs(H))-1) = payload_data(i:i+length(abs(H))-1).*exp(-1j*AngleVector)./abs(H);
        k = k+1;
        i = i + conf.N;
        end
        end
end



if strcmp(conf.plotfigure,'true')
    figure(5);
    plot(real(payload_data),imag(payload_data),'bo');
    title('filtered constelation');
    grid on
    xlabel('Re');
    ylabel('Im');

end 

rxbits = demapper(payload_data);


end

 


