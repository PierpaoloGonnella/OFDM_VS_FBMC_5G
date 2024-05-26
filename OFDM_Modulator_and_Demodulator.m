% OFDM_Modulator_and_Demodulator.md

function out = OFDM_Modulator_P(mappedSignal, NumFFT, CPlen, NumSC)
    % Calculate block size
    block_size = size(mappedSignal, 1) * size(mappedSignal, 2) / NumSC;
    
    % Serial to Parallel conversion
    S2P = reshape(mappedSignal, block_size, NumSC);
    Sub_carrier = S2P; % Serial to Parallel conversion
    
    % IFFT Block
    ifft_Subcarrier = [];
    for j = 1:NumSC
        if(NumFFT > 0)
            ifft_Subcarrier(:,j) = ifft(S2P(:,j), NumFFT);
        end
        if(NumFFT == 0)
            ifft_Subcarrier(:,j) = ifft(S2P(:,j));
        end
    end
    
    % Add Cyclic Prefix
    cp_start = block_size - CPlen;
    cyclic_prefix = [];
    Append_prefix = [];
    for i = 1:NumSC
        for j = 1:CPlen
            cyclic_prefix(j, i) = ifft_Subcarrier(j + cp_start, i);
        end
        Append_prefix(:, i) = vertcat(cyclic_prefix(:, i), ifft_Subcarrier(:, i));
    end
    
    % Conversion to serial stream for transmission
    [rows_Append_prefix, cols_Append_prefix] = size(Append_prefix);
    len_ofdm_data = rows_Append_prefix * cols_Append_prefix;
    
    % OFDM signal to be transmitted
    out = reshape(Append_prefix, 1, len_ofdm_data);
end

function out = OFDM_Demodulator_P(receivedSignal, NumFFT, CPlen, NumSC)
    % Calculate block size
    rxBlock_size = size(receivedSignal, 1) * size(receivedSignal, 2) / NumSC;
    
    % Remove cyclic Prefix
    rx_Subcarrier = reshape(receivedSignal, rxBlock_size, NumSC);
    rx_Subcarrier(1:CPlen, :) = [];
    
    % FFT Block
    fft_data = [];
    for w = 1:NumSC
        if(NumFFT > 0)
            fft_data(:, w) = fft(rx_Subcarrier(:, w), NumFFT);
        end
        if(NumFFT == 0)
            fft_data(:, w) = fft(rx_Subcarrier(:, w));
        end
    end
    
    % Conversion to serial and demodulation
    size_serial_data = [size(fft_data, 2) * size(fft_data, 1), 1];
    out = reshape(fft_data, size_serial_data(1), 1);
end
