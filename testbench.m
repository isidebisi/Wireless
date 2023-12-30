%% Wireless receiver testbench
%%RUN AUDIOTRANS SEVERAL TIMES TO FIND THE BEST CONFIGURATION
clear, clc, close all
conf = conf();


nbcarrier = [128 256 512 1024];
nbcarrierDefault = 256;
berNbCarrier = [];
txlenNbCarrier = [];

spacing = [10 5 2 1];
spacingDefault = 1;
berSpacing = [];
txlenSpacing = [];

f_train = [1 3 5 10 20 40];
f_trainDefault = 10;
berF_train = [];
txlenF_train = [];

cp_len = [256 128 32 1];
cp_lenDefault = 32;
berCp_Len = [];
txlenCp_len = [];

for ii = nbcarrier
    [a b] = audiotrans(ii,spacingDefault,f_trainDefault, cp_lenDefault, conf);
    berNbCarrier = [berNbCarrier a];
    txlenNbCarrier = [txlenNbCarrier b];
    if a <= min(berNbCarrier)
        nbcarrierDefault = a;
    end
end

    
for ii = spacing
    [a b] = audiotrans(nbcarrierDefault,ii,f_trainDefault, cp_lenDefault, conf);
    berSpacing = [berSpacing a];
    txlenSpacing = [txlenSpacing b];
    if a <= min(berSpacing)
        spacingDefault = a;
    end
end

for ii = f_train
    [a b] = audiotrans(nbcarrierDefault,spacingDefault,ii, cp_lenDefault, conf);
    berF_train = [berF_train a];
    txlenF_train = [txlenF_train b];
    if a <= min(berF_train)
        f_trainDefault = a;
    end
end

for ii = cp_len
    [a b] = audiotrans(nbcarrierDefault,spacingDefault,f_trainDefault, ii, conf);
    berCp_Len = [berCp_Len a];
    txlenCp_len = [txlenCp_len b];
    if a <= min(berCp_Len)
        cp_lenDefault = a;
    end
end

%% Plotting
figure(1)
subplot(2,2,1)
plot(nbcarrier, berNbCarrier)
title('BER vs. Number of Carriers')
xlabel('Number of Carriers')
ylabel('BER')

subplot(2,2,2)
plot(spacing, berSpacing)
title('BER vs. Spacing')
xlabel('Spacing')
ylabel('BER')

subplot(2,2,3)
plot(f_train, berF_train)
title('BER vs. Training Frequency')
xlabel('Training Frequency')
ylabel('BER')

subplot(2,2,4)
plot(cp_len, berCp_Len)
title('BER vs. Cyclic Prefix Length')
xlabel('Cyclic Prefix Length')
ylabel('BER')

