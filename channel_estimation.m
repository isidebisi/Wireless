function corrected_data = channel_estimation(rx_data, conf)

    payload_data = zeros(conf.nbcarrier * conf.OFDM_symbols, 1);
    symbol_index = 1; % Index for symbol processing

    for k = 1:size(rx_data, 2) % Iterate through columns

        switch conf.estimation_type    
            case 'none'
                if mod(k, conf.f_train + 1) == 1
                    continue;
                else
                    payload_data(symbol_index:symbol_index + conf.nbcarrier - 1) = osfft(rx_data(:,k), conf.os_factor);
                    symbol_index = symbol_index + conf.nbcarrier;
                end

            case 'viterbi'
                if mod(k, conf.f_train + 1) == 1
                    Y = osfft(rx_data(:,k), conf.os_factor);
                    H = Y ./ conf.train_seq;
                    AngleVector = angle(Y ./ conf.train_seq);
                    PreviousAngle = angle(H);
                else
                    payload_data(symbol_index:symbol_index + conf.nbcarrier - 1) = osfft(rx_data(:,k), conf.os_factor);
                    Segment = payload_data(symbol_index:symbol_index + conf.nbcarrier - 1);
                    for j = 1:conf.nbcarrier
                       A = pi / 2 * (-1:4);
                       deltaTheta = 1 / 4 * angle(Segment(j)^4) + A;
                       [~, ind] = min(abs(deltaTheta - PreviousAngle(j)));
                       AngleVector(j) = deltaTheta(ind);
                       AngleVector(j) = mod(0.01 * AngleVector(j) + 0.99 * PreviousAngle(j), 2 * pi);
                       PreviousAngle(j) = AngleVector(j);
                    end
                    payload_data(symbol_index:symbol_index + length(abs(H)) - 1) = payload_data(symbol_index:symbol_index + length(abs(H)) - 1) .* exp(-1j * AngleVector) ./ abs(H);
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








