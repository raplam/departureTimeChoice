function congestion = generateBottleneck(S)
% Very simple function creating a constant-capacity bottleneck.
% Last modified by Raphael Lamotte, on October 24, 2018.
congestion.mechanism='bottleneck';
congestion.S=S; % capacity
end

