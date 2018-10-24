This MATLAB code simulates the day-to-day adjustment of departure time decisions within a population. The demand is considered inelastic (its size is given) and users only have to schedule one trip (typically the morning or evening commute). When choosing their departure time, users attempt to maximize their utility, which depends both on their departure and arrival times. Congestion arises when many users have similar schedule preferences.

The code allows for different schedule preference functions, different adjustment mechanisms, and most importantly, for two congestion mechanisms: a constant capacity bottleneck and a Macroscopic Fundamental Diagram (MFD). It is also possible to describe the population either as a continuum, or with discrete agents.

When the population is described by a continuum, the set of possible departure has to be finite and specified in advance (for instance, users can depart at every minute within a 4 h time period). When the population is discrete however, users can depart at any time.

This repository includes several files:
“main.m” allows launching any single simulation. It includes options that can be modified by users.
The 3 files starting by “script” contains the scripts used to generate the figures of the 3rd chapter and conclusion of my PhD thesis, entitled “Congestion and departure time choice equilibrium in urban road networks”.
The other folders contain various pieces of code that are called by the scripts above.

Note: The combination “discrete agents” & “bottleneck” is not supported yet. I do not know whether it will ever be, because it does not seem very useful to me (the bottleneck model is more naturally implemented with a continuum of users).
