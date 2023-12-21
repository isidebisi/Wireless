function conf= conf()

conf.audiosystem = 'bypass';          
conf.plotfigure = 'true';
conf.estimation_type = 'viterbi';


% OFDM
conf.SNR = 30;
conf.f_s     = 48000;   % sampling rate  
conf.nbcarrier = 256;
conf.f_sym   = 5/2;     % symbol rate (K/2 where K is a divisor of 3x5^3)
conf.nframes =1;       % number of frames to transmit
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 8000;
conf.npreamble  = 100;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.f_train = 10; % Training symbol placed after every N data symbols
conf.rolloff_factor = 0.22;
conf.os_factor = conf.f_s/(conf.f_sym*conf.nbcarrier);
conf.filterlength = 5*conf.os_factor;
conf.figurecount = 1;
conf.enlarged_bandwidth = floor((conf.nbcarrier+1)/2)*conf.f_sym*1.1;


conf.cp_len = floor(conf.nbcarrier*0.2); 
conf.train_seq = -2*(randi([0 1],conf.nbcarrier,1)) + 1;

end

