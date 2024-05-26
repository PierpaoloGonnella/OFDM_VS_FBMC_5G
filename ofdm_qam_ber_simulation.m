% OFDM with QAM Modulation Bit Error Rate Simulation

% Parameters
numSC = 512;           % Number of OFDM subcarriers
cpLen = 64;            % OFDM cyclic prefix length
maxBitErrors = 1000;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted
M = 16;                % Modulation order

% Initialize plot
figure

% Initialize variables
k = 1;  

% OFDM Modulator and Demodulator setup
ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);

% AWGN Channel setup
channel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

% Error Rate Calculator setup
errorRate = comm.ErrorRate('ResetInputPort',true);

% Get OFDM Modulator dimensions
ofdmDims = info(ofdmMod);
numDC = ofdmDims.DataInputSize(1);
frameSize = [k*numDC 1];

% Calculate Eb/No vectors for SNR
EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(log2(M));

% Initialize BER vector and error statistics
berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);  
    
% Iterate over SNR values
for m = 1:length(EbNoVec)
    snr = snrVec(m);
    
    % Simulation loop
    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        % Generate binary data
        dataIn = randi([0,1],frameSize);
       
        % QAM modulation
        qamTx = qammod(dataIn,M,'bin');
        
        % OFDM modulation
        txSig = ofdmMod(qamTx);
        
        % Calculate Tx signal power
        powerDB = 10*log10(var(txSig));
        
        % Calculate the noise variance
        noiseVar = 10.^(0.1*(powerDB-snr));
        
        % Pass the signal through a noisy channel
        rxSig = channel(txSig,noiseVar);
        
        % OFDM demodulation
        qamRx = ofdmDemod(rxSig);
        
        % QAM demodulation
        dataOut = qamdemod(qamRx,M,'bin'); 
        
        % Collect error statistics
        errorStats = errorRate(dataIn,dataOut,0);
    end
    
    % Save BER data
    berVec(m,:) = errorStats;
    
    % Reset the error rate calculator
    errorStats = errorRate(dataIn,dataOut,1);
end

% Theoretical BER calculation
berTheory = berawgn(EbNoVec,'qam',M);

% Plot results
semilogy(EbNoVec,berVec(:,1),'^') % Simulated BER
hold on
semilogy(EbNoVec,berTheory)       % Theoretical BER
hold on
grid on

% Reset and release modulator and demodulator objects
reset(ofdmMod);
reset(ofdmDemod);
release(ofdmMod);
release(ofdmDemod);
