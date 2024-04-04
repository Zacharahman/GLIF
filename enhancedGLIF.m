function [V, spikeTimes] = enhancedGLIF(n, simDuration, dt)
    % Basic neuron parameters
    E_L = -65; % Resting potential (mV)
    V_th = -50; % Threshold potential (mV)
    V_reset = -65; % Reset potential (mV)
    tau_m = 20; % Membrane time constant (ms)
    R_m = 1; % Membrane resistance (MΩ)
    I_ext = 1.5; % External current (nA)
    noise_sigma = 2; % Standard deviation for noise
    
    % Initialize state variables
    t = 0:dt:simDuration; % Time vector
    V = E_L * ones(n, length(t)); % Membrane potential matrix
    spikeTimes = cell(n, 1); % Cell array to record spike times
    
    % Define initial network connectivity
    W = randn(n) * 0.2; % Synaptic weights matrix
    W(logical(eye(size(W)))) = 0; % Eliminate self-connections
    
    % Define synaptic plasticity parameters
    alpha_increase = 0.005; % Increase in synaptic weight for each spike
    
    % Simulation loop
    for i = 2:length(t)
        for neuron = 1:n
            % Calculate synaptic current based on current synaptic weights
            I_syn = sum(W(neuron, :) .* (V(:, i-1) > V_th)') * R_m;
            
            % Integrate membrane potential with noise
            I_noise = noise_sigma * randn;
            dV = ((E_L - V(neuron, i-1)) + R_m * (I_ext + I_syn) + I_noise) * dt / tau_m;
            V(neuron, i) = V(neuron, i-1) + dV;
            
            % Check for spike event
            if V(neuron, i) >= V_th
                V(neuron, i) = V_reset; % Reset after spike
                spikeTimes{neuron} = [spikeTimes{neuron}, t(i)]; % Record spike time
                
                % Implement synaptic plasticity (enhancement)
                W(:, neuron) = W(:, neuron) + alpha_increase; % Strengthen synapses after spiking
            end
        end
    end
    
    % Visualization of membrane potentials for each neuron
    for neuron = 1:n
        figure; plot(t, V(neuron, :));
        title(sprintf('Neuron %d Membrane Potential', neuron));
        xlabel('Time (ms)'); ylabel('Membrane potential (mV)');
    end
end
