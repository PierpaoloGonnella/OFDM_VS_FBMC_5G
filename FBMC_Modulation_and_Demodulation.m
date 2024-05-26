% FBMC_Modulation_and_Demodulation.m

% FBMC Modulator: Prototype Filter, Upsampling, Filtering, and IFFT

function out = FBMC_Modulator_P(mappedSignal, numFFT, numGuards, K, numSymbols)
% mappedSignal: Input mapped symbols
% numFFT: Number of FFT points
% numGuards: Number of guard subcarriers
% K: Number of subcarriers per symbol
% numSymbols: Number of symbols

% Prototype Filter
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end

% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% Transmit-end processing
% Initialize arrays
L = numFFT - 2 * numGuards;  % Number of complex symbols per OFDM symbol
KF = K * numFFT;
KL = K * L;
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);
sumFBMCSpec = zeros(KF * 2, 1);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2 * KF, 1));

% Loop over symbols
for symIdx = 1:numSymbols
    % OQAM Modulator: alternate real and imaginary parts
    if rem(symIdx,2)==1     % Odd symbols
        dataSubCar(1:2:L) = real(mappedSignal);
        dataSubCar(2:2:L) = 1i * imag(mappedSignal);
    else                    % Even symbols
        dataSubCar(1:2:L) = 1i * imag(mappedSignal);
        dataSubCar(2:2:L) = real(mappedSignal);
    end
    
    % Upsample by K, pad with guards, and filter with the prototype filter
    dataSubCarUp(1:K:end) = dataSubCar;
    dataBitsUpPad = [zeros(numGuards * K, 1); dataSubCarUp; zeros(numGuards * K, 1)];
    X1 = filter(Hk, 1, dataBitsUpPad);
    % Remove 1/2 filter length delay
    X = [X1(K:end); zeros(K - 1, 1)];

    % Compute IFFT of length KF for the transmitted symbol
    txSymb = fftshift(ifft(X));

    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2 + 1:end); complex(zeros(numFFT/2, 1))];
    symBuf(KF + (1:KF)) = symBuf(KF + (1:KF)) + txSymb;

    % Compute power spectral density (PSD)
    currSym = complex(symBuf(1:KF));
    [specFBMC, fFBMC] = periodogram(currSym, hann(KF, 'periodic'), KF * 2, 1);
    sumFBMCSpec = sumFBMCSpec + specFBMC;

    % Store transmitted signals for all symbols
    txSigAll(:, symIdx) = currSym;
end
     
% Plot power spectral density
sumFBMCSpec = sumFBMCSpec / mean(sumFBMCSpec(1 + K + 2 * numGuards * K:end - 2 * numGuards * K - K));
plot(fFBMC - 0.5, 10 * log10(sumFBMCSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['FBMC, K = ' num2str(K) ' overlapped symbols'])
set(gcf, 'Position', figposition([15 50 30 30]));

out = currSym;
end


% FBMC Demodulator: FFT, Filtering, and OQAM Processing

function out = FBMC_Demodulator_P(receivedSignal, numGuards, K, numSymbols)
% receivedSignal: Received FBMC signal
% numGuards: Number of guard subcarriers
% K: Number of subcarriers per symbol
% numSymbols: Number of symbols

% Prototype Filter Coefficients (Assumed to be known)
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end

% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% Process symbol-wise
for symIdx = 1:numSymbols
    rxSig = receivedSignal(:, symIdx);

    % Perform FFT
    rxf = fft(fftshift(rxSig));

    % Matched filtering with prototype filter
    rxfmf = filter(Hk, 1, rxf);
    % Remove K-1 delay elements
    rxfmf = [rxfmf(K:end); zeros(K - 1, 1)];
    % Remove guards
    rxfmfg = rxfmf(numGuards * K + 1:end - numGuards * K);

    % OQAM post-processing
    %  Downsample by 2K, extract real and imaginary parts
    if rem(symIdx, 2)
        % Imaginary part is K samples after real one
        r1 = real(rxfmfg(1:2 * K:end));
        r2 = imag(rxfmfg(K + 1:2 * K:end));
        rcomb = complex(r1, r2);
    else
        % Real part is K samples after imaginary one
        r1 = imag(rxfmfg(1:2 * K:end));
        r2 = real(rxfmfg(K + 1:2 * K:end));
        rcomb = complex(r2, r1);
    end
    %  Normalize by the upsampling factor
    rcomb = (1 / K) * rcomb;

end
    
out = rcomb;
end
