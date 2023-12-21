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
diff = conf.nbits - conf.requiredBits
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

train_ifft = osifft(conf.train_seq,conf.os_factor);
train_ifft = [train_ifft(end-conf.cp_len+1:end);train_ifft];

%% Prepare the symbols in OFDM symbol format
symbolMatrix = reshape(tx_symbols,[conf.nbcarrier,conf.OFDM_symbols]).';

TimeMatrix = zeros(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train), conf.f_s/conf.spacing+conf.cp_len);
k = 1;

% add Training symbols every so often for continuous frame estimation
for i=1:(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))
    if mod(i,conf.f_train+1) == 0
        TimeMatrix(i,:) =train_ifft;
    else
    TimeSignalTemp = osifft(symbolMatrix(k,:),conf.os_factor);
    TimeMatrix(i,1:conf.cp_len) = TimeSignalTemp(end-conf.cp_len+1:end);
    TimeMatrix(i,conf.cp_len+1:end) = TimeSignalTemp;
    k = k+1;
    end
end

tx_ifft = reshape(TimeMatrix.',(conf.OFDM_symbols+ floor(conf.OFDM_symbols/conf.f_train))*(conf.f_s/conf.spacing+conf.cp_len),1);


%% normalizing the signals
preamble_norm = preamble_upsample_filtered/sqrt( mean(real(preamble_upsample_filtered).^2) + mean(imag(preamble_upsample_filtered).^2));
train_ifft_norm = train_ifft/sqrt( mean(real(train_ifft).^2) + mean(imag(train_ifft).^2));
tx_ifft_norm = tx_ifft/sqrt( mean(real(tx_ifft).^2) + mean(imag(tx_ifft).^2));

txsignalTemp = [preamble_norm; train_ifft_norm; tx_ifft_norm];

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



