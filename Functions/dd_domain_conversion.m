function yDD = dd_domain_conversion(y_t, N, M, Ts, F0, dt)
    yDD = zeros(M, N);
    for n = 0:N-1
        f = n * F0;
        exp_t = exp(1j * 2 * pi * f * (0:dt:(length(y_t)-1) * dt));
        y_conv = fft_convolve(y_t, exp_t, dt);
        
        % Sample at delay taps
        t_taps = (0:M-1) * Ts;
        [~, loc] = min(abs(t_taps' - (0:dt:(length(y_t)-1)*dt)), [], 2);
        yDD(:,n+1) = y_conv(loc);
    end
    yDD = yDD(:);
end
