% Zach Abdulrahman | Generalized Leaky Integrate Fire Model

function [spikeTimes, V] = GLIF(n, simDuration, dt)
    % Parameters
    E_L = -65; % Resting potential (mV)
    V_th = -50; % Threshold potential (mV)
    V_spike = 30; % Spike potential (mV)
    V_reset = -65; % Reset potential (mV)
    tau_m = 20; % Membrane time constant (ms)
    R_m = 1; % Membrane resistance (MΩ)
    I_ext = 1.5; % External current (nA)
    noise_sigma = 2; % Standard deviation for noise
    
    % Initialize state variables
    t = 0:dt:simDuration; % Time vector
    V = E_L * ones(n, length(t)); % Membrane potential matrix
    spikeTimes = cell(n, 1); % Cell array to record spike times
    
    % Random synaptic weight matrix (W_ij: effect of neuron j on neuron i)
    W = randn(n) * 0.2; % Random weights for demonstration
    for i = 1:n
        W(i, i) = 0; % No self-connections
    end
    
    % Simulation loop
    for i = 2:length(t)
        for neuron = 1:n
            if V(neuron, i-1) == V_spike
                V(neuron, i) = V_reset; % Reset after spike
            else
                I_syn = sum(W(neuron, :) .* (V(:, i-1) > V_th)') * R_m; % Synaptic current
                I_noise = noise_sigma * randn; % Gaussian noise
                dV = ((E_L - V(neuron, i-1)) + R_m * I_ext + I_syn + I_noise) * dt / tau_m;
                V(neuron, i) = V(neuron, i-1) + dV;
                
                % Spike event
                if V(neuron, i) >= V_th
                    V(neuron, i) = V_spike;
                    spikeTimes{neuron} = [spikeTimes{neuron}, t(i)];
                end
            end
        end
    end
    
  % Start of the plotting section for each neuron
    for neuron = 1:n
        figure; % Create a new figure for each neuron
        plot(t, V(neuron, :)); % Plot the membrane potential over time
        hold on; % Hold on to plot spikes on the same figure
        spikes = spikeTimes{neuron}; % Retrieve spike times for this neuron
        for j = 1:length(spikes) % Loop through each spike time
            plot([spikes(j), spikes(j)], [E_L, V_spike], 'k', 'LineWidth', 2); % Mark the spike
        end
        hold off; % Release the figure for normal operations
        title(sprintf('Neuron %d', neuron)); % Title with neuron number
        xlabel('Time (ms)'); % Label for the x-axis
        ylabel('Membrane potential (mV)'); % Label for the y-axis
    end
end


%n = 5; % Number of neurons
%simDuration = 1000; % Simulation duration in milliseconds
%dt = 1; % Time step in milliseconds
%[spikeTimes, V] = GLIF(n, simDuration, dt);

