% OFDM_Spectrum_Analysis.md

% Parameters
F = 1000;
Ts = 1 / F;
M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 4;             % Number of OFDM subcarriers
cpLen = 4;             % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 512;      % Maximum number of bits transmitted
numGuards = 212;

% Number of FFT points and complex symbols per OFDM symbol
numFFT = 4096;
L = numFFT - 2 * numGuards;

% QPSK Modulator and Demodulator setup
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

sumOFDMSpec = 0;

% Loop over symbols
for symIdx = 1:maxNumBits / 4
    % Generate random binary data
    inpData2 = randi([0 1], 2 * L, 1);
    
    % Modulate data with QPSK
    modData = qpskMod(inpData2);

    % OFDM symbol construction
    symOFDM = [zeros(numGuards, 1); modData; zeros(numGuards, 1)];
    
    % IFFT and normalization
    ifftOut = sqrt(numFFT) .* ifft(ifftshift(symOFDM));

    % Power spectral density (PSD) calculation
    [specOFDM, fOFDM] = periodogram(ifftOut, rectwin(length(ifftOut)), ...
        numFFT * 2, 1, 'centered');
    sumOFDMSpec = sumOFDMSpec + specOFDM;
end

% Plot power spectral density (PSD) over all subcarriers
sumOFDMSpec = sumOFDMSpec / mean(sumOFDMSpec(1 + 2 * numGuards:end - 2 * numGuards));
figure;
plot(fOFDM, 10 * log10(sumOFDMSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['OFDM, numFFT = ' num2str(numFFT)])
set(gcf, 'Position', figposition([46 50 30 30]));
