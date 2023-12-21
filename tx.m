function [txsignal conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

%% Map txbits into Symbols using QPSK

if conf.nbits ~= conf.requiredBits
diff = conf.nbits - conf.requiredBits;
added = randi([0 1],diff,1);
txbits = [txbits;added]; % padding with random bits to form a vector of length a multiple of 512


end


tx_symbols = mapGray(txbits);

%% Generate Preamble into BPSK


preamble = -2 * preambleGenerate(conf.npreamble)+1;
preamble_upsample = upsample(preamble, conf.os_factor);
preamble_upsample_filtered = matched_filter(preamble_upsample, conf);
preamble_upsample_filtered = preamble_upsample_filtered(conf.filterlength+1:end-conf.filterlength);

%% Prepare the training data

Training_time = osifft(conf.train_seq,conf.os_factor);
Training_Vector = [Training_time(end-conf.cp_len+1:end);Training_time];

%% Prepare the symbols in OFDM symbol format
symbolMatrix = reshape(tx_symbols,[conf.nbcarrier,conf.OFDM_symbols]).';

TimeMatrix = zeros(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train), conf.f_s/conf.f_sym+conf.cp_len);
k = 1;
for i=1:(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))
    if mod(i,conf.f_train+1) == 0
        TimeMatrix(i,:) =Training_Vector;
    else
    TimeSignalTemp = osifft(symbolMatrix(k,:),conf.os_factor);
    TimeMatrix(i,1:conf.cp_len) = TimeSignalTemp(end-conf.cp_len+1:end);
    TimeMatrix(i,conf.cp_len+1:end) = TimeSignalTemp;
    k = k+1;
    end
end

TimeVector = reshape(TimeMatrix.',(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))*(conf.f_s/conf.f_sym+conf.cp_len),1);


%% normalizing the signals
preamble_upsample_filtered = preamble_upsample_filtered/sqrt( mean(real(preamble_upsample_filtered).^2) + mean(imag(preamble_upsample_filtered).^2));
Training_Vector = Training_Vector/sqrt( mean(real(Training_Vector).^2) + mean(imag(Training_Vector).^2));
TimeVector = TimeVector/sqrt( mean(real(TimeVector).^2) + mean(imag(TimeVector).^2));

txsignalTemp = [preamble_upsample_filtered;Training_Vector; TimeVector];

if strcmp(conf.plotfigure,'true')
    fftSignal = abs(fftshift(fft(txsignalTemp)));
    N = length(fftSignal); 
    f = (-N/2:N/2-1) * (conf.f_s/N);

    figure(6);
    plot(f,fftSignal);
    title('Spectrum of transmitted signal before modulation');
    xlabel('frequency (Hz)');
    ylabel('Amplitude');
    
   

end


%% Modulate
txsignal = zeros(length(txsignalTemp),1);
t = 1/conf.f_s :1/conf.f_s:length(txsignal)/conf.f_s;
txsignal = real(txsignalTemp) .* cos(2*pi*conf.f_c*t).' - imag(txsignalTemp).*sin(2*pi*conf.f_c*t).';



