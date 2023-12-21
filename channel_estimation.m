function corrected_data = channel_estimation(RX_Time_Matrix, conf)
    TrainingData = conf.train_seq;

    payload_data = zeros(conf.nbcarrier * conf.OFDM_symbols, 1);
    k = 1; % OFDM symbol index
    i = 1; % QPSK symbol index

    while i < length(payload_data)
        switch conf.estimationtype    
            case 'none'
                if mod(k, conf.f_train + 1) == 1
                    k = k + 1;
                else
                    payload_data(i:i + conf.nbcarrier - 1) = osfft(RX_Time_Matrix(conf.cp_len + 1:end, k), conf.os_factor);
                    k = k + 1;
                    i = i + conf.nbcarrier;
                end

            case 'viterbi'
                if mod(k, conf.f_train + 1) == 1
                    Y = osfft(RX_Time_Matrix(conf.cp_len + 1:end, k), conf.os_factor);
                    H = Y ./ TrainingData;
                    AngleVector = angle(H);
                    PreviousAngle = angle(H);
                    k = k + 1;
                else
                    payload_data(i:i + conf.nbcarrier - 1) = osfft(RX_Time_Matrix(conf.cp_len + 1:end, k), conf.os_factor);
                    Segment = payload_data(i:i + conf.nbcarrier - 1);
                    for j = 1:conf.nbcarrier
                       A = pi / 2 * (-1:4);
                       deltaTheta = 1 / 4 * angle(Segment(j)^4) + A;
                       [~, ind] = min(abs(deltaTheta - PreviousAngle(j)));
                       AngleVector(j) = deltaTheta(ind);
                       AngleVector(j) = mod(0.01 * AngleVector(j) + 0.99 * PreviousAngle(j), 2 * pi);
                       PreviousAngle(j) = AngleVector(j);
                    end
                    payload_data(i:i + length(abs(H)) - 1) = payload_data(i:i + length(abs(H)) - 1) .* exp(-1j * AngleVector) ./ abs(H);
                    k = k + 1;
                    i = i + conf.nbcarrier;
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





