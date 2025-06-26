function x_hat = ODDM_equalizer(yDD,HDD,N,M,L_range,Es,N0,S,N_iters,RwDD)

% Find all possible Lambda_n matrices and Theta_n matrices - Looks right
possible_Lambda_m = zeros(N*M,N);
possible_Theta_m = zeros(N*M,N);
for m = 0:M-1

    % Find second indices
    indices2 = m*N+1:(m+1)*N;

    % Find new Lambda_m and Theta_m
    Lambda_m = zeros(N);
    Theta_m = zeros(N);
    for l = 0:M-1
        if ismember(l-m,L_range)
            % Find first indices for the first H block
            indices11 = l*N+1:(l+1)*N;

            % Second first H block
            selected_H_block1 = HDD(indices11,indices2);

            % Add to Lambda_m
            Lambda_m_add = selected_H_block1' * selected_H_block1;
            Lambda_m = Lambda_m + Lambda_m_add;

            for l_prime = 0:M-1
                if ismember(l_prime-m,L_range)
                    % Find first indices for the second H block
                    indices12 = l_prime*N+1:(l_prime+1)*N;

                    % Second second H block and RwDD block
                    selected_R_block = RwDD(indices11,indices12);
                    selected_H_block2 = HDD(indices12,indices2);

                    % Add to Lambda_m
                    Theta_m_add = selected_H_block1' * selected_R_block * selected_H_block2;
                    Theta_m = Theta_m + Theta_m_add;
                end
            end
        end
    end

    % Add Lambda_n and Theta_n to stack
    possible_Lambda_m(indices2,:) = Lambda_m;
    possible_Theta_m(indices2,:) = Theta_m;
end

% possible_Theta_m = possible_Lambda_m;
1;

% Predefine variables and start iterator equalizer
x_hat = zeros(N*M,1);
flag_detector = true;
iters = 0;
while iters < N_iters && flag_detector
    iters = iters + 1;

    % Sweep through all M blocks of yDD (size Nx1)
    for l = 0:M-1

        % Find first indices for current block
        indices1 = l*N+1:(l+1)*N;

        % Create y_tilde for current l
        y_l = yDD(indices1);

        % Create gamma_n with both for loops
        gamma_m = zeros(N,1);

        for m = 0:M-1
            if ismember(l-m,L_range)

                % Calculate ISI from chosen x_hat_l
                ISI = zeros(N,1);
                for m_prime = 0:M-1
                    if ismember(l-m_prime,L_range) && m_prime ~= m

                        % Find second indices for current block
                        indices2 = m_prime*N+1:(m_prime+1)*N;

                        % Select current blocks and make ISI to add
                        selected_H_block = HDD(indices1,indices2);
                        selected_x_block = x_hat(indices2);
                        ISI_add = selected_H_block * selected_x_block;

                        % Add to ISI
                        ISI = ISI + ISI_add;
                    end
                end

                % Remove noise from y_hat_l
                y_hat_l = y_l - ISI;

                % Find second indices for current block
                indices2 = m*N+1:(m+1)*N;

                % Add to Gamma_n
                selected_H_block = HDD(indices1,indices2);
                gamma_m = gamma_m + selected_H_block' * y_hat_l;
            end
        end

        % Select Lambda_n from possible matrices
        Lambda_m = possible_Lambda_m(indices1,:);
        Theta_m = possible_Theta_m(indices1,:);

        % Create MMSE matrix for iterative solver
        W_n = Lambda_m' / (Lambda_m*Lambda_m' + (N0/Es)*Theta_m);

        % Create x_hat for current block and push to stack
        x_hat_n = W_n * gamma_m;

        % Hard detection for block
        dist1 = abs(x_hat_n.' - S).^2;
        [~,min_index1] = min(dist1);
        x_hat_n = S(min_index1);

        % Map hard encoded x_hat_n to x_hat
        x_hat(indices1) = x_hat_n;
    end

    % Check if duplicate result is found, break if true
    if iters > 1
        if last_x_hat == x_hat
            flag_detector = false;
        end
    end

    % Save last x_hat
    last_x_hat = x_hat;
end
