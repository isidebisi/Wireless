function conf= conf()

conf.audiosystem = 'bypass';          
conf.plotfigure = 'false';
conf.estimation_type = 'viterbi';


% OFDM
conf.f_s     = 48000;   % sampling rate  
conf.nbcarrier = 128;
conf.spacing   = 5;    
conf.nframes =1;       % number of frames to transmit
conf.f_c     = 8000;
conf.npreamble  = 100;
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.f_train = 3; % Training symbol placed after every N data symbols
conf.rolloff_factor = 0.22;
conf.os_factor = conf.f_s/(conf.spacing*conf.nbcarrier);
conf.filterlength = 5*conf.os_factor;
conf.figurecount = 1;
conf.enlarged_bandwidth = floor((conf.nbcarrier+1)/2)*conf.spacing*1.2;


conf.cp_len = 30; 
conf.train_seq = -2*(preambleGenerate(conf.nbcarrier)) + 1;

end

