function RwDD = gen_RwDD(N,M,T,shape,alpha,Q,res)

% Add redundancy for rectangular pulses
if shape == "rect"
    Q = 1;
end

% Define file and folder names
file_name = sprintf('DD_table_N%d_M%d_T%d_%s_a%.2f_Q%d_%dres',N,M,T,shape,alpha,Q,res);
folder_name = 'DD Noise Covariance Tables';
full_path = folder_name + "/" + file_name + ".mat";

% Check if the folder exists, and create it if it does not
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% Define ambiguity table
try

    % Load file if it exists
    load(full_path)

    1;

catch

    fprintf("File not found for delay-Doppler noise covariance matrix...\nGenerating...\n")

    % Set parameters
    Ts = T / M;
    F0 = 1 / (N*T);
    tau = (M:-1:-M)*Ts;

    % Render all ambiguity values
    RwDD = zeros(N*M);
    update_vals = floor((N*N)*(linspace(.01,.99,34)));
    fprintf("Progress:                        |\n")
    count = 0;
    ambig_vec = zeros(1,length(tau));
    for n1 = 0:N-1
        for n2 = 0:N-1

            % Update progress bar
            count = count + 1;
            if ismember(count,update_vals)
                fprintf("x")
            end

            % Set frequency instances
            f1 = n1*F0;
            f2 = n2*F0;

            % Render vector of values for this block
            for i = 1:length(tau)
                ambig_vec(i) = Apq(tau(i),f1-f2,N,M,T,shape,alpha,Q,res) * exp(1j*2*pi*f2*tau(i));
            end

            % Generate toeplitz block
            block = zeros(M);
            for i = 0:M-1
                this_row = circshift(ambig_vec,i+M+1);
                block(i+1,:) = this_row(1:M);
            end

            % Assign block to matrix
            RwDD((n1*M+1):((n1+1)*M),(n2*M+1):((n2+1)*M)) = block;

        end
    end
    fprintf("\nComplete!\n\n")

    % Save to be used later
    save(full_path,"RwDD")

end