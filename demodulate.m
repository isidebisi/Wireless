function demodulated_signal = demodulate(modulated_signal, conf)
% demodulator
%   
%   modulated_signal : input signal
%   conf : global configuration variable 
%
%   demodulated_signal : output demodulated signal

    t = 0:1/conf.f_s:(length(modulated_signal)-1)/conf.f_s;

    % carrier wave
    carrier_wave = exp(-1i*2*pi*(conf.f_c*t'));
    demodulated_signal = modulated_signal .* carrier_wave;
end


