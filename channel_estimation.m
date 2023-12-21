function corrected_data = channel_estimation(rx_data, conf)

    payload_data = zeros(conf.nbcarrier * conf.ofdm_symbols, 1);
    symbol_index = 1; % Index for symbol processing

    for k = 1:size(rx_data, 2) % Iterate through columns

        switch conf.estimation_type    
            case 'none'
                % remove training sequence between symbols when there is
                % one
                if mod(k, conf.f_train + 1) == 1
                    continue;
                else
                    % training sequence removal and frequency domain
                    % conversion of symbols
                    payload_data(symbol_index:symbol_index + conf.nbcarrier - 1) = osfft(rx_data(:,k), conf.os_factor);
                    symbol_index = symbol_index + conf.nbcarrier;
                end

            case 'viterbi'
                if mod(k, conf.f_train + 1) == 1
                    % frequency domain conversion 
                    Y = osfft(rx_data(:,k), conf.os_factor);

                    H = Y ./ conf.train_seq;
                    vect = angle(Y ./ conf.train_seq);
                    prev_angle = angle(H);
                else
                    % frequency domain conversion of symbol
                    payload_data(symbol_index:symbol_index + conf.nbcarrier - 1) = osfft(rx_data(:,k), conf.os_factor);
                    seg = payload_data(symbol_index:symbol_index + conf.nbcarrier - 1);

                    for j = 1:conf.nbcarrier
                   
                       % difference between two  successive phase offset is
                       % a gaussian random variable
                       deltaTheta = 1 / 4 * angle(seg(j)^4) + pi / 2 * (-1:4) ;
                       
                       [~, ind] = min(abs(deltaTheta - prev_angle(j)));
                       vect(j) = deltaTheta(ind);
                       vect(j) = mod(0.01 * vect(j) + 0.99 * prev_angle(j), 2 * pi);
                       prev_angle(j) = vect(j);
                    end
                    payload_data(symbol_index:symbol_index + length(abs(H)) - 1) = payload_data(symbol_index:symbol_index + length(abs(H)) - 1) .* exp(-1j * vect) ./ abs(H);
                    symbol_index = symbol_index + conf.nbcarrier;
                end
        end
    end

    if strcmp(conf.plotfigure, 'true')
        figure(5);
        plot(real(payload_data), imag(payload_data), 'bo');
        title('filtered constellation');
        grid on
        xlabel('Re');
        ylabel('Im');
    end

    corrected_data = payload_data;
end








