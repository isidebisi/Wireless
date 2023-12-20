% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework by Maroun Wakim and Jan Clevorn
%
%   4 operating modes for audiosystem:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal
%   - 'awgn'   : Applies simultaed awgn and phase noise
%
%   2 type of data to be sent:
%   - 'random' : random generated bitstream
%   - 'image'  : sends an image
%
%   3 type of channel estimation
%   - 'block'  : uses block estimation
%   - 'viterbi': uses viterbi-viterbi algorithm to estimate
%   - 'none'   : no channel estimation
%
%   Plotfigure to view plots fot analysis:
%   - 'true'   : view all plots
%   - 'false'  : view no plots
%
% Configuration Values
% % % % %

conf.audiosystem = 'bypass';    
conf.datatype = 'image';       
conf.plotfigure = 'true';
conf.estimationtype = 'viterbi';




% OFDM
conf.SNR = 30;
conf.sigmaDeltaTheta = 0.004;
conf.f_s     = 48000;   % sampling rate  
conf.N = 256;
conf.f_sym   = 5/2;     % symbol rate (K/2 where K is a divisor of 3x5^3)
conf.nframes =1;       % number of frames to transmit
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 8000;
conf.npreamble  = 100;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.repeatTrainingFrequency = 10; % Training symbol placed after every N data symbols
conf.rolloff_factor = 0.22;
conf.os_factor = conf.f_s/(conf.f_sym*conf.N);
conf.filterlength = 5*conf.os_factor;
conf.figurecount = 1;
conf.enlarged_bandwidth = floor((conf.N+1)/2)*conf.f_sym*1.1;


if mod(conf.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

for k=1:conf.nframes

    if strcmp(conf.datatype,'random')
        % Generate random data
        conf.requiredBits = 2000;
        conf.nbits   = conf.requiredBits + 2*conf.N - mod(conf.requiredBits,2*conf.N);    % number of bits 
        txbits = randi([0 1],conf.requiredBits,1);

    elseif strcmp(conf.datatype,'image')
        % load and show image
        YourImage = imread('girl.jpeg');
        if strcmp(conf.plotfigure,'true')
            figure(1)
            hold on
            imshow(YourImage); title('Transmitted Original Image');
        end 
        
        orig_class = class(YourImage);
        orig_size = size(YourImage);
        txbits = reshape((dec2bin(typecast(YourImage(:), 'uint8'), 8) - '0').', 1, []).';
        conf.requiredBits = length(txbits);
        conf.nbits   = conf.requiredBits + 2*conf.N - mod(conf.requiredBits,2*conf.N);    % number of bits 
    end
    
    conf.OFDM_symbols = conf.nbits/conf.modulation_order/conf.N;
    conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

    % Transmit Function
    [txsignal conf] = tx(txbits,conf,k);
    
    % % % % % % % % % % % %
    % Begin
    % Audio Transmission
    %
    
    % normalize values
    peakvalue       = max(abs(txsignal));
    normtxsignal    = txsignal / (peakvalue + 0.3);
    
    % create vector for transmission
    rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
    rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
    
    if strcmp(conf.plotfigure,'true')
        figure(2);
        plot(rawtxsignal(:,1));
        title('Transmited Tx signal');
        xlabel('Sample');
        ylabel('Amplitude');
    end 

    %wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
    audiowrite('out.wav',rawtxsignal,conf.f_s)  
    
    % Platform native audio mode 
    if strcmp(conf.audiosystem,'native')
        
        % Windows WAV mode 
        if ispc()
            disp('Windows WAV');
            wavplay(rawtxsignal,conf.f_s,'async');
            disp('Recording in Progress');
            rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
            disp('Recording complete')
            rxsignal = rawrxsignal(1:end,1);

        % ALSA WAV mode 
        elseif isunix()
            disp('Linux ALSA');
            cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
            system(cmd); 
            disp('Recording in Progress');
            system('aplay  out.wav')
            pause(2);
            disp('Recording complete')
            rawrxsignal = audioread('in.wav');
            rxsignal    = rawrxsignal(1:end,1);
        end
        
    % MATLAB audio mode
    elseif strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;

    elseif strcmp(conf.audiosystem,'awgn')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;
    end
    

    if strcmp(conf.plotfigure,'true')
        figure(conf.figurecount + 1);
        plot(rxsignal);
        title('Received Signal')

    end 

    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    % Receive Function
    [rxbits, conf]       = rx(rxsignal,conf);
    rxbits = rxbits(1:end - (conf.nbits - conf.requiredBits),:);
    res.rxnbits(k)      = length(rxbits);  
    res.biterrors(k)    = sum(rxbits ~= txbits);
    
    

end


if strcmp(conf.datatype,'image')
    reconstructed = reshape(typecast(uint8(bin2dec(char(reshape(rxbits, 8, [])+'0').')), orig_class), orig_size);

    if strcmp(conf.plotfigure,'true')
        figure(3);
        imshow(reconstructed); title('Recieved Reconstructured Image')

    end 
end


per = sum(res.biterrors > 0) / conf.nframes;
ber = sum(res.biterrors) / sum(res.rxnbits);

disp(['Estimation Type: ', num2str(conf.estimationtype)]);
disp(['Packet Error Rate (PER): ', num2str(per)]);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

