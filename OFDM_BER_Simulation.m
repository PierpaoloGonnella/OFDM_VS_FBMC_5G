% OFDM_BER_Simulation.m

% OFDM Bit Error Rate (BER) Simulation

% Parameters
M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 512;           % Number of OFDM subcarriers
cpLen = 96;            % OFDM cyclic prefix length
maxBitErrors = 10^37;  % Maximum number of bit errors
maxNumBits = 2048;     % Maximum number of bits transmitted
Nifft = 128;           % IFFT size

% QPSK Modulator and Demodulator setup
qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

% AWGN Channel setup
channel = comm.AWGNChannel('NoiseMethod','Variance','VarianceSource','Input port');

% Error Rate Calculator setup
errorRate = comm.ErrorRate('ResetInputPort',true);

numDC = 4; % Number of data carriers

frameSize = [k*numDC 1]; % Frame size for serial-to-parallel conversion

numPAPR = 1; % Number of PAPR measurements

ofdmTheory = [];
ofdmMeasurement = [];

% Eb/No and SNR vectors setup
EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/Nifft);

% BER vectors initialization
berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);

% BER Simulation
for m = 1:length(EbNoVec)
    snr = snrVec(m);
    
    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        % Generate binary data
        dataIn = randi([0,1],maxNumBits,1);
        
        % Apply QPSK modulation
        qpskTx = qpskMod(dataIn);
        
        % OFDM Modulation
        ofdmMod = OFDM_Modulator_P(qpskTx, 0, cpLen, Nifft);
        
        % Calculate Tx signal power
        powerDB = 10*log10(var(ofdmMod));
        
        % Calculate the noise variance
        noiseVar = 10.^(0.1*(powerDB-snr)); 
        
        % Pass the signal through a noisy channel
        afterChannel = channel(ofdmMod, noiseVar);
        
        % OFDM Demodulation
        recvdSerialData = OFDM_Demodulator_P(afterChannel, 0, cpLen, numSC);
        
        % Apply QPSK Demodulation
        qpskRx = qpskDemod(recvdSerialData);
     
        % Collect error statistics
        errorStats = errorRate(dataIn, qpskRx, 0);
    end
    
    % Save BER data
    berVec(m,:) = errorStats;   
    
    % Reset the error rate calculator
    errorStats = errorRate(dataIn, qpskRx, 1);   
    
    % PAPR calculation
    lenOFDMData = size(ofdmMod, 1) * size(ofdmMod, 2);
    r2 = abs(ofdmMod);
    modx = reshape(r2, 1, lenOFDMData);
    sigma2 = mean(modx.^2, 1);
    PAPR(m) = max(modx.^2, 1) / sigma2;
    
    % Append data for CCDF calculation
    ofdmMeasurement(:,numPAPR) = reshape(afterChannel, size(afterChannel, 2), 1);
    ofdmTheory(:,numPAPR) = reshape(ofdmMod, size(ofdmMod, 2) * size(ofdmMod, 1), 1);
end

% Theoretical BER calculation
berTheory = berawgn(EbNoVec,'psk',M,'nondiff');

% CCDF calculation
[CCDFy, CCDFx, PAPR] = ccdf(ofdmMeasurement(:,1));
plot(ccdf); hold on;
[CCDFy, CCDFx, PAPR] = ccdf(ofdmTheory(:,1));

xlabel('PAPR'); 
ylabel('CCDF');
title('OFDM Simulation');
legend('QPSK', 'Location', 'Best')
grid on;

hold off;

% OFDM Modulator function
function out = OFDM_Modulator_P(mappedSignal, NumFFT, CPlen, NumSC)
    block_size = size(mappedSignal, 1) * size(mappedSignal, 2) / NumSC;
    S2P = reshape(mappedSignal, block_size, NumSC);
    
    Sub_carrier = S2P; % Serial to Parallel conversion
    
    ifft_Subcarrier = [];
    for j = 1:NumSC
        if(NumFFT > 0)
            ifft_Subcarrier(:,j) = ifft(S2P(:,j), NumFFT);
        end
        if(NumFFT == 0)
            ifft_Subcarrier(:,j) = ifft(S2P(:,j));
        end
    end
    
    cp_start = block_size - CPlen;
    cyclic_prefix = [];
    Append_prefix = [];
    
    for i = 1:NumSC
        for j = 1:CPlen
            cyclic_prefix(j, i) = ifft_Subcarrier(j+cp_start, i);
        end
        Append_prefix(:, i) = vertcat(cyclic_prefix(:, i), ifft_Subcarrier(:, i));
    end
    
    len_ofdm_data = size(Append_prefix, 1) * size(Append_prefix, 2);
    out = reshape(Append_prefix, 1, len_ofdm_data);
end

% OFDM Demodulator function
function out = OFDM_Demodulator_P(receivedSignal, NumFFT, CPlen, NumSC)
    rxBlock_size = size(receivedSignal, 1
    ) * size(receivedSignal, 2) / NumSC;
    rx_Subcarrier = reshape(receivedSignal, rxBlock_size, NumSC);
    
    rx_Subcarrier(1:CPlen, :) = [];
    
    fft_data = [];
    for w = 1:NumSC
        if(NumFFT > 0)
            fft_data(:, w) = fft(rx_Subcarrier(:, w), NumFFT);
        end
        if(NumFFT == 0)
            fft_data(:, w) = fft(rx_Subcarrier(:, w));
        end
    end

    size_serial_data = [size(fft_data, 2) * size(fft_data, 1), 1];
    out = reshape(fft_data, size_serial_data(1), 1);
end
