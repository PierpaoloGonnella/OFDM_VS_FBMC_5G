function run_ofdm_modulator_tests()
% Unit tests for the OFDM modulator defined in OFDM_Modulator_and_Demodulator.m.
% Runs in both MATLAB and GNU Octave (plain assert-based, no toolboxes needed).
%
%   addpath(repoRoot); addpath(fullfile(repoRoot, 'tests'));
%   run_ofdm_modulator_tests
%
% The file's first function (OFDM_Modulator_P) is reachable through the file
% name, so the modulator is called as OFDM_Modulator_and_Demodulator(...).

tests = {@test_output_length_and_shape, ...
         @test_cyclic_prefix_copies_block_tail, ...
         @test_fft_of_body_recovers_input, ...
         @test_zero_nfft_matches_explicit_nfft, ...
         @test_parseval_energy, ...
         @test_linearity, ...
         @test_impulse_gives_constant_symbol};
failed = 0;
for i = 1:numel(tests)
    name = func2str(tests{i});
    try
        tests{i}();
        fprintf('PASS  %s\n', name);
    catch err
        failed = failed + 1;
        fprintf('FAIL  %s: %s\n', name, err.message);
    end
end
fprintf('%d/%d tests passed\n', numel(tests) - failed, numel(tests));
if failed > 0
    error('ofdm:tests', '%d test(s) failed', failed);
end
end

function out = mod_ofdm(x, nfft, cp, nsc)
out = OFDM_Modulator_and_Demodulator(x, nfft, cp, nsc);
end

function x = rand_qpsk(n, seed)
rand('seed', seed); %#ok<RAND>
b = rand(n, 2) > 0.5;
x = ((2*b(:,1) - 1) + 1i*(2*b(:,2) - 1)) / sqrt(2);
end

function test_output_length_and_shape()
N = 16; cp = 4; nsc = 8;
y = mod_ofdm(rand_qpsk(N*nsc, 1), N, cp, nsc);
assert(isequal(size(y), [1, (N + cp) * nsc]));
end

function test_cyclic_prefix_copies_block_tail()
N = 16; cp = 4; nsc = 8;
y = mod_ofdm(rand_qpsk(N*nsc, 2), N, cp, nsc);
blocks = reshape(y, N + cp, nsc);
for k = 1:nsc
    assert(max(abs(blocks(1:cp, k) - blocks(end-cp+1:end, k))) < 1e-12);
end
end

function test_fft_of_body_recovers_input()
N = 32; cp = 8; nsc = 4;
x = rand_qpsk(N*nsc, 3);
y = mod_ofdm(x, N, cp, nsc);
blocks = reshape(y, N + cp, nsc);
rx = fft(blocks(cp+1:end, :));
assert(max(abs(rx(:) - x(:))) < 1e-10);
end

function test_zero_nfft_matches_explicit_nfft()
N = 16; cp = 2; nsc = 4;
x = rand_qpsk(N*nsc, 4);
assert(max(abs(mod_ofdm(x, 0, cp, nsc) - mod_ofdm(x, N, cp, nsc))) < 1e-12);
end

function test_parseval_energy()
% IFFT scales energy by 1/N per block (body only, CP excluded).
N = 64; cp = 16; nsc = 2;
x = rand_qpsk(N*nsc, 5);
blocks = reshape(mod_ofdm(x, N, cp, nsc), N + cp, nsc);
body = blocks(cp+1:end, :);
assert(abs(sum(abs(body(:)).^2) - sum(abs(x).^2) / N) < 1e-10);
end

function test_linearity()
N = 16; cp = 4; nsc = 4;
a = rand_qpsk(N*nsc, 6); b = rand_qpsk(N*nsc, 7);
lhs = mod_ofdm(2*a - 3i*b, N, cp, nsc);
rhs = 2*mod_ofdm(a, N, cp, nsc) - 3i*mod_ofdm(b, N, cp, nsc);
assert(max(abs(lhs - rhs)) < 1e-12);
end

function test_impulse_gives_constant_symbol()
% A single active subcarrier at DC yields a flat time-domain block of 1/N.
N = 8; cp = 2; nsc = 1;
x = zeros(N, 1); x(1) = 1;
y = mod_ofdm(x, N, cp, nsc);
assert(max(abs(y - 1/N)) < 1e-12);
end
