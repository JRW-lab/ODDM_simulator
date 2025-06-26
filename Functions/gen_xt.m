function [x_t,g_t] = gen_xt(xDD, N, M, T, shape, alpha, Q, res)

Ts = T / M;
F0 = 1 / (N*T);
dt = Ts / res;
xDD_block = reshape(xDD, M, N);
t_range = 0:dt:(N*T + Q*Ts);

% Generate sample filter to get normalization factor
a_t_sample = gen_pulse(0:dt:Q*Ts,shape,Ts,Q,alpha);
norm_val = sqrt(sum(abs(a_t_sample).^2) * dt) * sqrt(N);

n1 = (0:(N-1)).';
g_t = zeros(M,length(t_range));
x_t = zeros(1,length(t_range));
for m = 0:M-1

    % Set up TX pulse train g(t - m*Ts)
    t_range_at = (t_range - m*Ts) - n1*T;
    a_t = gen_pulse(t_range_at,shape,Ts,Q,alpha);
    g_t(m+1,:) = sum(a_t, 1) / norm_val;


    for n = 0:N-1

        exp_t = exp(1j * 2 * pi * n * F0 * (t_range - m*Ts));

        x_t = x_t + xDD_block(m+1,n+1) .* g_t(m+1,:) .* exp_t;

    end

end