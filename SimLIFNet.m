function [spk, NetParams, V] = SimLIFNet(W,varargin)
simTime = 100;
tstep = 0.01;
initialConditions = [];
refractoryTime = 0;
offsetCurrents = 0;
synapticDensity = 4;
forcingFunctions = {};
noiseAmplitude = 0;
displayProgress = 1;
plotResults = 1;

% Parse varargin
if nargin > 1
    if isstruct(varargin{1})
        NetParams = varargin{1};
    else
        paramNames = varargin(1:2:end);
        paramValues = varargin(2:2:end);
        for iParam = 1:length(paramNames)
            switch paramNames{iParam}
                case 'simTime'
                    simTime = paramValues{iParam};
                case 'tstep'
                    tstep = paramValues{iParam};
                case 'initialConditions'
                    initialConditions = paramValues{iParam};
                case 'refractoryTime'
                    refractoryTime = paramValues{iParam};
                case 'offsetCurrents'
                    offsetCurrents = paramValues{iParam};
                case 'synapticDensity'
                    synapticDensity = paramValues{iParam};
                case 'forcingFunctions'
                    forcingFunctions = paramValues{iParam};
                case 'noiseAmplitude'
                    noiseAmplitude = paramValues{iParam};
                case 'displayProgress'
                    displayProgress = paramValues{iParam};
                case 'plotResults'
                    plotResults = paramValues{iParam};
                otherwise
                    error('Unknown parameter: %s', paramNames{iParam});
            end
        end
    end
end

% Simulation setup
nNeurons = length(W);
nSteps = round(simTime / tstep);
spk = cell(nNeurons, 1);
V = zeros(nNeurons, nSteps);
if isempty(initialConditions)
    V(:, 1) = rand(nNeurons, 1);
else
    V(:, 1) = initialConditions;
end
t = 0;
NetParams = struct('W', W, 'simTime', simTime, 'tstep', tstep, ...
    'initialConditions', initialConditions, 'refractoryTime', refractoryTime, ...
    'offsetCurrents', offsetCurrents, 'synapticDensity', synapticDensity, ...
    'forcingFunctions', forcingFunctions, 'noiseAmplitude', noiseAmplitude, ...
    'displayProgress', displayProgress, 'plotResults', plotResults);

% Main simulation loop
for iStep = 2:nSteps
    % Compute synaptic inputs
    Isynapse = zeros(nNeurons, 1);
    for j = 1:nNeurons
        if ~isempty(spk{j})
            Isynapse(j) = sum(synapticDensity(j) * exp(-synapticDensity(j) * (t - spk{j})) .* W(:, j));
        end
    end
    
    % Compute forcing inputs
    for i = 1:length(forcingFunctions)
        idx = forcingFunctions{i, 2};
        forcingFunction = forcingFunctions{i, 1};
        Isynapse(idx) = Isynapse(idx) + forcingFunction(t);
    end
    
    % Compute noise
    noise = noiseAmplitude * randn(nNeurons, 1);
    
    % Update membrane potential
    V(:, iStep) = V(:, iStep - 1) + tstep * (-V(:, iStep - 1) + offsetCurrents + Isynapse + noise);
    
    % Check for spikes
    for j = 1:nNeurons
        if V(j, iStep) >= 1
            spk{j} = [spk{j}, t];
            V(j, iStep) = 0;
        end
    end
    
    % Apply refractory period
    for j = 1:nNeurons
        if ~isempty(spk{j}) && t - spk{j}(end) < refractoryTime
            V(j, iStep) = 0;
        end
    end
    
    % Update time
    t = t + tstep;
    
    % Display progress
    if displayProgress
        fprintf('Simulation Progress: %.2f%%\r', 100 * (iStep / nSteps));
    end
end

% Plot results
if plotResults
    figure;
    subplot(2, 1, 1);
    plot((1:nSteps) * tstep, V);
    title('Membrane Voltage');
    xlabel('Time');
    ylabel('Voltage');
    subplot(2, 1, 2);
    for j = 1:nNeurons
        plot(spk{j}, j * ones(size(spk{j})), '.');
        hold on;
    end
    title('Raster Plot');
    xlabel('Time');
    ylabel('Neuron');
    hold off;
end
