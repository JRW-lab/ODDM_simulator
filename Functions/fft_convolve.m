function y = fft_convolve(x, h, dt)
    N = length(x) + length(h) - 1;
    Y = ifft(fft(x, N) .* fft(h, N)) * dt;
    y = Y(1:length(x)); % Truncate to original length
end
