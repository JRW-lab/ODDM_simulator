function x_shifted = fractional_delay(x, tau, dt)
    % Create a sinc-based fractional delay filter
    delay_samples = tau / dt;
    h = sinc((0:20) - delay_samples);
    h = h / sum(h); % Normalize filter energy
    x_shifted = conv(x, h, 'same'); % Apply filter
end
