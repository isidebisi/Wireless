% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Project
% Charlotte Heibig and Ismael Frei
%
clear, clc, close all

conf = conf();

if mod(conf.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

for k=1:conf.nframes


    % load and show image
    image = imread('resizeLudovic.png');
    image = rgb2gray(image);
    if strcmp(conf.plotfigure,'true')
        figure(1)
        hold on
        imshow(image); title('Transmitted Original Image');
    end 
    

    imageSize = size(image);

    txbits = reshape((dec2bin(typecast(image(:), 'uint8'), 8) - '0').', 1, []).';
    conf.requiredBits = length(txbits);
    conf.nbits   = conf.requiredBits + 2*conf.N - mod(conf.requiredBits,2*conf.N);    % number of bits 

    
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



reconstructed = reshape(typecast(uint8(bin2dec(char(reshape(rxbits, 8, [])+'0').')), 'uint8'), imageSize);

if strcmp(conf.plotfigure,'true')
    figure(3);
    imshow(reconstructed); title('Recieved Reconstructured Image')

end


per = sum(res.biterrors > 0) / conf.nframes;
ber = sum(res.biterrors) / sum(res.rxnbits);


disp(['Packet Error Rate (PER): ', num2str(per)]);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

