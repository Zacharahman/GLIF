% Ask the user for the number of neurons for the GLIF simulation
n = input('Enter the number of neurons (n): ');

% Fixed simulation parameters
simDuration = 1000; % Fixed simulation duration in milliseconds
dt = 1; % Fixed time step in milliseconds

% Validate the input for the number of neurons
if isempty(n) || n <= 0
    error('Number of neurons must be a positive integer.');
end

% Assuming the GLIF function is defined in another file 'GLIF.m' and is accessible
[spikeTimes, V] = GLIF(n, simDuration, dt);

% Notify the user that the simulation is complete
disp('Simulation complete. Check the figures for results.');
