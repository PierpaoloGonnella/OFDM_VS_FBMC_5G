% Parameters for QPSK
M = 4;                  % Modulation alphabet for QPSK
numSC = 512;            % Number of OFDM subcarriers
cpLen = 64;             % OFDM cyclic prefix length
maxBitErrors = 1000;    % Maximum number of bit errors
maxNumBits = 1e7;       % Maximum number of bits transmitted

% Additional Parameters for 16-QAM
M16QAM = 16;            % Modulation alphabet for 16-QAM
EbNoVec16QAM = (0:10)'; % Eb/No Vector for 16-QAM
snrVec16QAM = EbNoVec16QAM + 10*log10(log2(M16QAM)); % SNR Vector for 16-QAM

% Initialize figure
figure

% QPSK Modulator and Demodulator
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

% OFDM Modulator and Demodulator for QPSK
ofdmMod = comm.OFDMModulator('FFTLength', numSC, 'CyclicPrefixLength', cpLen);
ofdmDemod = comm.OFDMDemodulator('FFTLength', numSC, 'CyclicPrefixLength', cpLen);

% AWGN Channel
channel = comm.AWGNChannel('NoiseMethod', 'Variance', ...
    'VarianceSource', 'Input port');

% Error Rate Calculator for QPSK
errorRateQPSK = comm.ErrorRate('ResetInputPort', true);

% Get OFDM dimensions
ofdmDims = info(ofdmMod);
numDC = ofdmDims.DataInputSize(1);
frameSize = [log2(M) * numDC 1];

% Eb/No Vector calculation for QPSK
EbNoVecQPSK = (0:10)';
snrVecQPSK = EbNoVecQPSK + 10*log10(log2(M)) + 10*log10(numDC/numSC);

% BER Vector initialization for QPSK
berVecQPSK = zeros(length(EbNoVecQPSK), 3);
errorStatsQPSK = zeros(1, 3);

% Simulation loop for different Eb/No values for QPSK
for m = 1:length(EbNoVecQPSK)
    snr = snrVecQPSK(m);
    
    while errorStatsQPSK(2) <= maxBitErrors && errorStatsQPSK(3) <= maxNumBits
        % Generate binary data
        dataIn = randi([0, 1], frameSize);
        
        % Apply modulation and demodulation for QPSK
        qpskTx = qpskMod(dataIn);
        txSig = ofdmMod(qpskTx);
        powerDB = 10*log10(var(txSig));
        noiseVar = 10^(0.1*(powerDB-snr));
        rxSig = channel(txSig, noiseVar);
        qpskRx = ofdmDemod(rxSig);
        dataOut = qpskDemod(qpskRx);
        
        % Collect error statistics for QPSK
        errorStatsQPSK = errorRateQPSK(dataIn, dataOut, 0);
    end
    
    % Save BER data for QPSK
    berVecQPSK(m, :) = errorStatsQPSK;
    errorStatsQPSK = errorRateQPSK(dataIn, dataOut, 1);
end

% Theoretical BER calculation for QPSK
berTheoryQPSK = berawgn(EbNoVecQPSK, 'qam', M);

% Plotting for QPSK
semilogy(EbNoVecQPSK, berTheoryQPSK, '^')
hold on
semilogy(EbNoVecQPSK, berTheoryQPSK)

% Reset and release OFDM objects for QPSK
reset(ofdmMod);
reset(ofdmDemod);
release(ofdmMod);
release(ofdmDemod);

% Error Rate Calculator for 16-QAM
errorRate16QAM = comm.ErrorRate('ResetInputPort', true);

% BER Vector initialization for 16-QAM
berVec16QAM = zeros(length(EbNoVec16QAM), 3);
errorStats16QAM = zeros(1, 3);

% Simulation loop for different Eb/No values for 16-QAM
for m = 1:length(EbNoVec16QAM)
    snr = snrVec16QAM(m);
    
    while errorStats16QAM(2) <= maxBitErrors && errorStats16QAM(3) <= maxNumBits
        % Generate binary data
        dataIn = randi([0, 1], frameSize);
        
        % Apply modulation and demodulation for 16-QAM
        qamTx = qammod(dataIn, M16QAM, 'bin');
        txSig = ofdmMod(qamTx);
        powerDB = 10*log10(var(txSig));
        noiseVar = 10^(0.1*(powerDB-snr));
        rxSig = channel(txSig, noiseVar);
        qamRx = ofdmDemod(rxSig);
        dataOut = qamdemod(qamRx, M16QAM, 'bin');
        
        % Collect error statistics for 16-QAM
        errorStats16QAM = errorRate16QAM(dataIn, dataOut, 0);
    end
    
    % Save BER data for 16-QAM
    berVec16QAM(m, :) = errorStats16QAM;
    errorStats16QAM = errorRate16QAM(dataIn, dataOut, 1);
end

% Theoretical BER calculation for 16-QAM
berTheory16QAM = berawgn(EbNoVec16QAM, 'qam', M16QAM);

% Plotting for 16-QAM
semilogy(EbNoVec16QAM, berTheory16QAM, 'o')
semilogy(EbNoVec16QAM, berTheory16QAM)
legend('Simulation QPSK', 'Theory QPSK', 'Simulation 16-QAM', ...
    'Theory 16-QAM', 'Location', 'Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off
