function filtered_signal = matched_filter(signal, conf)
h = rrc(conf.os_factor, conf.rolloff_factor, conf.filterlength);

filtered_signal = conv(h, signal,"full");