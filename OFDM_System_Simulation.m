% OFDM_System_Simulation.md

M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 4096*4;        % Number of OFDM subcarriers
cpLen = 8;             % OFDM cyclic prefix length
maxBitErrors = 10^37;  % Maximum number of bit errors
maxNumBits = 512;      % Maximum number of bits transmitted
Nifft = 128;

qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);
ccdf = comm.CCDF('PAPROutputPort', true, 'MaximumPowerLimit', 50);
errorRate = comm.ErrorRate('ResetInputPort', true);
numDC = 4;
frameSize = [k*numDC, 1];
iN = 256;
numPAPR = 1;
step = 0.1;
ofdmTeory = [];
ofdmMesurement = [];

while iN < numSC + 1
    EbNoVec = (0:10)';
    snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/iN);
    berVec = zeros(length(EbNoVec), 3);
    errorStats = zeros(1, 3);
    
    for m = 1:length(EbNoVec)
        snr = snrVec(m);
        while errorStats(2) <= maxBitErrors && errorStats(3) <= 4096
            dataIn = randi([0, 1], maxNumBits, 1);  
            dataInVector = reshape(dataIn, maxNumBits*maxNumBits, 1);
            qpskTx = qammod(dataInVector, 16, 'bin');
            ofdmMod = OFDM_Modulator_P(qpskTx, 0, cpLen, iN);
            if(size(qpskTx, 1) > 2)
                S2P = reshape(qpskTx, size(qpskTx, 1)/iN, iN);
            end
            Sub_carrier = [];
            for l = 1:iN
                Sub_carrier(:, l) = S2P(:, l);
            end
            ifft_Subcarrier = [];
            for j = 1:iN
                if(NumFFT > 0)
                    ifft_Subcarrier(:, j) = ifft(S2P(:, j), NumFFT);
                end
                if(NumFFT == 0)
                    ifft_Subcarrier(:, j) = ifft(S2P(:, j));
                end
            end 
            cp_start = block_size - cpLen;
            cyclic_prefix = [];
            Append_prefix = [];
            for i = 1:NumSC
                for j = 1:cpLen
                    cyclic_prefix(j, i) = ifft_Subcarrier(j + cp_start, i);
                end
                Append_prefix(:, i) = vertcat(cyclic_prefix(:, i), ifft_Subcarrier(:, i));
            end
            [rows_Append_prefix, cols_Append_prefix] = size(Append_prefix);
            len_ofdm_data = size(ifft_Subcarrier, 2) * size(ifft_Subcarrier, 1);
            out = reshape(Append_prefix, 1, len_ofdm_data);
            r2 = abs(out);
            modx = reshape(r2, 1, len_ofdm_data);
            sigma2 = mean(modx.^2, 1);
            txSig = ofdmMod;
            powerDB = 10*log10(var(txSig));
            noiseVar = 10.^(0.1*(powerDB-snr));
            if(maxNumBits > 64)
                channel = comm.AWGNChannel('NoiseMethod', 'Variance', 'EbNo', EbNoVec(1, m), 'VarianceSource', 'Input port');
            else
                channel = comm.AWGNChannel('NoiseMethod', 'Variance', 'VarianceSource', 'Input port');
            end
            after_channel = channel(ofdm_signal, noiseVar);
            rs1 = reshape(after_channel, size(after_channel, 2) * size(after_channel, 1), 1);
            rs2 = reshape(rs1, size(rs1, 1) / iN, iN);
            r2s = abs(after_channel);
            modxs = reshape(r2s, 1, size(after_channel, 2) * size(after_channel, 1));
            sigma2s = mean(modxs.^2, 1);
            recvd_serial_data = OFDM_Demodulator_P(after_channel, 0, cpLen, numSC);
            qpskRx = qamdemod(recvd_serial_data, 16, 'bin');
            d1 = reshape(dataIn, maxNumBits * maxNumBits, 1);
            errorStats = errorRate(d1, d1, 0);
        end
        berVec(m, :) = errorStats;
        errorStats = errorRate(d1, d1, 1);
        PAPR(iN, m) = max(modx.^2, 1) / sigma2;
        PAPRs(iN, m) = max(modxs.^2, 1) / sigma2s;
    end
    x = PAPR(iN, :);
    ofdmMesurement(:, numPAPR) = reshape(after_channel, size(after_channel, 2), 1);
    ofdmTeory(:, numPAPR) = reshape(ofdm_signal, size(ofdm_signal, 2), 1);
    numPAPR = numPAPR + 1;
    iN = iN * 4;
    berTheory = berawgn(EbNoVec, 'psk', M, 'nondiff');
end

[CCDFy, CCDFx, PAPR] = ccdf([ofdmMesurement(:, 2), ofdmMesurement(:, 3)]);
plot(ccdf); hold on;
[CCDFy, CCDFx, PAPR] = ccdf([ofdmTeory(:, 2), ofdmTeory(:, 3)]);
plot(ccdf); hold on;
xlabel('PAPR'); 
ylabel('CCDF');
title(' ');
legend('256 NFFT', '1024 NFFT');
grid on;
hold off;

function out = OFDM_Modulator_P(mappedSignal, NumFFT, CPlen, NumSC)
    block_size = size(mappedSignal, 1) * size(mappedSignal, 2) / NumSC;
    S2P = reshape(mappedSignal, block_size, NumSC);
    Sub_carrier = S2P;
    ifft_Subcarrier = [];
    for j = 1:NumSC
        if(NumFFT > 0)
            ifft_Subcarrier(:, j) = ifft(S2P(:, j), NumFFT);
        end
        if(NumFFT == 0)
            ifft_Subcarrier(:, j) = ifft(S2P(:, j));
        end
    end
    cp_start = block_size - CPlen;
   
    cyclic_prefix = [];
    Append_prefix = [];
    for i = 1:NumSC
        for j = 1:CPlen
            cyclic_prefix(j, i) = ifft_Subcarrier(j + cp_start, i);
        end
        Append_prefix(:, i) = vertcat(cyclic_prefix(:, i), ifft_Subcarrier(:, i));
    end
    [rows_Append_prefix, cols_Append_prefix] = size(Append_prefix);
    len_ofdm_data = rows_Append_prefix * cols_Append_prefix;
    out = reshape(Append_prefix, 1, len_ofdm_data);
end

function out = OFDM_Demodulator_P(receivedSignal, NumFFT, CPlen, NumSC)
    rxBlock_size = size(receivedSignal, 1) * size(receivedSignal, 2) / NumSC;
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
