%% Equalizer DD RX signal using iterative solver
clc; clear;
addpath(fullfile(pwd, 'Functions'));

% Simulation Settings
res = 10;
num_frames = 10;
EbN0_dB_vec = 3;

% Signal data
M_ary = 4;
Es = 1;
Eb = Es / log2(M_ary);

% Define inputs
N_iters = 10;
N = 16;
M = 64;
T = 1 / 15000;
Fc = 4e9;
v = 500;
% shape = "rect";
% shape = "sinc";
shape = "rrc";
alpha = 0.4;
Q = 8;

% Define parameters
syms_per_f = M*N;
Ts = T / M;
F0 = 1 / (N*T);

% Data setup
if M_ary == 2
    % Define bit order
    bit_order = [0;1];

    % Define alphabet
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    % Define bit order
    bit_order = [0,0;0,1;1,0;1,1];

    % Define alphabet
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

% Render ambiguity table
[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,v,shape,alpha,Q,res);

% Generate delay-Doppler noise covariance and it's "half" version
RwDD = eye(N*M);
RwDD_half = eye(N*M);

% Loop through simulator
BER = zeros(length(EbN0_dB_vec),1);
SER = zeros(length(EbN0_dB_vec),1);
FER = zeros(length(EbN0_dB_vec),1);
for j = 1:length(EbN0_dB_vec)

    % Define SNR and noise variance
    EbN0_dB = EbN0_dB_vec(j);
    N0 = Eb / (10^(EbN0_dB / 10));

    % Reset bit errors for each SNR
    bit_errors = zeros(num_frames,syms_per_f*log2(M_ary));
    sym_errors = zeros(num_frames,syms_per_f);
    frm_errors = zeros(num_frames,1);
    for i = 1:num_frames

        fprintf("%.2f%% (%.2f%%)     ||     Simulating SNR %d/%d, frame %d/%d\n", ...
            (j/length(EbN0_dB_vec)) * 100, (i/num_frames) * 100, j, length(EbN0_dB_vec), i, num_frames);

        % Generate data
        [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);

        % Generate H matrix and channel information
        [HDD,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(T,N,M,Fc,v,Q,Ambig_Table);

        % Generate noise
        zDD = sqrt(N0/2) * RwDD_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));

        % Construct received signal - Discrete
        yDD = HDD * xDD + zDD;

        % Equalize received signal
        x_hat = OTFS_pulse_equalizer_AWGN(yDD,HDD,N,M,L2,-L1,Es,N0,S,N_iters);

        % Hard detection for final x_hat
        dist = abs(x_hat.' - S).^2;
        [~,min_index] = min(dist);
        RX_sym = min_index.' - 1;

        % Convert final RX_sym to RX_bit
        RX_bit = bit_order(RX_sym+1,:);

        % Error calculation
        diff_bit = TX_bit ~= RX_bit;
        diff_sym = TX_sym ~= RX_sym;
        bit_errors(i,:) = diff_bit(:).';
        sym_errors(i,:) = diff_sym(:).';
        if sum(sym_errors(i,:)) > 0
            frm_errors(i) = 1;
        end

    end

    % Calculate BER, SER and FER
    BER(j) = sum(bit_errors,"all") / (num_frames*N*M*log2(M_ary));
    SER(j) = sum(sym_errors,"all") / (num_frames*N*M);
    FER(j) = sum(frm_errors,"all") / (num_frames);

    % % Print results
    % fprintf("BER: " + BER(j) + "\n")
    % fprintf("SER: " + SER(j) + "\n")
    % fprintf("FER: " + FER(j) + "\n")

end

semilogy(EbN0_dB_vec,BER)
grid on
ylim([1e-6 1e-1])
xlim([3 18])
xticks(4:2:18)