
%   Max Mason
%   RangeMaster_WithInternalResistance_VA_Mason.m

%   state space simulation of a Range Master with implementation of voltage
%   supply internal resistance.

%   guitar pedals sounding better when running out of battery is well
%   documented. 
%   When a 9V battery gets low on charge this is seen as an
%   increase in internal resistance rather than a direct drop in open
%   circuit voltage. 
%   In this script the internal resistance of the battery
%   can be entered which is used in tandem with other known values and
%   equations 3.19 in "2019 Mac Porter - Virtual Analog Modelling of Guitar
%   Effects Circuits" to find I1 (the current through the battery) in terms
%   of the battery voltage. This can then be substituted into V = Vin - I1*R 
%   which is then solved implicitly to give the equation found in line 154.

%   current settings apply to audio input, if sinosoidal input required
%   stability may become an issue and Newtons settings iteration maximums
%   may need to be increased. Current sinosoidal input options match the
%   results in Porter's dissertation.

clear all
close all

%%%%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   flags
plotting = true;                        % plotting on/off
input = "audio";                        % input == "audio" / "sin" / "delta"
sound = true;                           % play output sounds on/off
SaveAudio = false;                      % save output to file on/off

%   input
filename = "InputSample_Mason.wav";     % file to be read if input == audio
drive = 0.5;                            % multiply audio sample to get more distortion
SR = 88.2e3;                            % sample rate [Hz]
f0 = 1000;                              % frequency [Hz]
amp = 0.3;                              % amplitude of input
dur = 0.1;                              % duration of sim [s]

%   Physical parameters

VDC = -9;                               % nominal open circuit voltage of the battery
Rs = 5000;                              % internal resistance of the battery (0: fully charged -> 5000: basically dead)

r = 470e3;                              % resistance of R[Ohms]
r2 = 68e3;
r3 = 10e3;
r4 = 3.9e3;
r5 = 1e6;

c = 4.7e-9;                             % capacitance of C [F]
c2 = 47e-6;
c3 = 10e-9;
Is = 2.029e-6;                          % diode saturation current (A)
Vt = 25.83e-3;                          % diode thermal voltage (V)
Ni = 1.752;                             % diode ideality factor
Bf = 108.1;
Br = 8.602;

% Newton's method settings
tol = 1e-9;                             % convergence tolerance
max_iter = 10;                          % max number of iterations
limlim = 10;                            % limit of damping iterations

%%%%% input creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if input == "sin" 
    
    numSamples = round(SR*dur);             % number of samples in sim
    tvec = dur*(0:numSamples-1)'/numSamples;% time vector [s]
    u = amp*sin(2*pi*f0*tvec);              % input vector
    NT = numSamples;                        % number of samples
    k = 1/SR;                               % time step
    
end

if input == "audio"
    
    [u,SR] = audioread(filename);           % read in audio
    u = u*drive;                            % scale audio
    k = 1/SR;                               % time step
    [NT,unused] = size(u);                  % number of samples
    tvec = (0:NT-1)*k;                      % time vector (for plotting)
    
end

if input == "delta"
    
    numSamples = round(SR*dur);             % number of samples in sim
    tvec = dur*(0:numSamples-1)'/numSamples;% time vector [s]
    u = [1;zeros(numSamples-1,1)];          % input vector
    NT = numSamples;                        % number of samples
    k = 1/SR;                               % time step
    
end

%%%%% Derived parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   scheme parameters
A = [(-r-r2)/(c*r*r2),0,0;0,-1/(r4*c2),0;0,0,-1/(c3*(r3+r5))];
B = [(r+r2)/(c*r*r2),-1/(r*c);0,0;0,1/(c3*(r3+r5))];
C = [-1/c,0;-1/c2,-1/c2;0,r3/(c3*(r3+r5))];
D = [1,1,0;1,0,-r3/(r3+r5)];
E = [-1,0;-1,r5/(r3+r5)];
F = [0,0;0,r3*r5/(r3+r5)];
L = [0,0,-r5/(r3+r5)];
M = [0,r5/(r3+r5)];
N = [0,r3*r5/(r3+r5)];

[Q, unused] = size(A);                  % size of identity matrix

%   derived scheme parameters
Ha = inv((2/k)*eye(Q) - A);
Hb = ((2/k)*eye(Q) + A);
K = D*Ha*C + F;

%   non-linear function
fv = @(v) [Is*((1/Bf)*(exp(v(1)/Vt) - 1) + (1/Br)*(exp(v(2)/Vt) - 1)); Is*((exp(v(1)/Vt) - 1) - ((1+Br)/Br)*(exp(v(2)/Vt) - 1))];

%   constituents of jacobian
g1_dv1 = @(v) (Is/(Bf*Vt))*exp(v(1)/Vt);
g1_dv2 = @(v) (Is/(Br*Vt))*exp(v(2)/Vt);
g2_dv1 = @(v) (Is/Vt)*exp(v(1)/Vt);
g2_dv2 = @(v) (Is/Vt)*((1+Br)/Br)*exp(v(2)/Vt);

%   initial values
x = [1.026426585775843;-0.904172291035517;-6.700314922196817];
x1 = [1.026426585775843;-0.904172291035517;-6.700314922196817];
v = [0.122254294740326;-5.673888336420975];
v1 = [0.122254294740326;-5.673888336420975];
y = zeros(NT,1);
I = [0.001870541203148;0.229968507780318]*0.001;
I1 = [0.001870541203148;0.229968507780318]*0.001;

u = [u,ones(NT,1)*VDC];                 % attach battery voltage to input

%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = [2:NT]                              % for each sample
    
num = 0;                                    % reset iteration counter
diff = 10;                                  % reset difference above tol

%   find battery output voltage from nominal voltage 
u(n,2) = (u(n,2) - u(n,2)/(r3/r + 1) - (Rs/r)*(x(1) - u(n,1)) + (v(2)*Rs/r)/(r3/r + 1)) / (1 - 1/(r3/r + 1) + Rs/r);

%   input to non-linearity
p = D*Ha*(Hb*x1 + B*u(n-1,:)' + C*I1) + (D*Ha*B + E)*u(n,:)';


while diff > tol && num < max_iter                          % while solution is outside of tolerance
    
    gv = K*fv(v)+p-v;                                       % equation to be solved
    J = (K*[g1_dv1(v),g1_dv2(v);g2_dv1(v),g2_dv2(v)]-1);    % jacobian
    DV = (J\gv);                                            % delta V
      
    v = v1 - DV;                                            % new guess
     
    lim = 1;                                                % reset damping iteration counter
    
    %   Damping loop, while next guess is bigger and while iteration limit not reached
    while max(abs(K*fv(v) + p - v)) > max(abs(K*fv(v1) + p - v1)) && lim < limlim
        
    v = v1 - (2^-lim)*DV;                                       % damped guess
      
    lim = lim + 1;                                              % step counter
    
    end
    
    diff = max(abs(v - v1),[],'all');                       % find difference in guesses
    num = num + 1;                                          % step iteration counter
    v1 = v;                                                 % pass current guesss to next step
    
end

I = fv(v);                                              % calculate current from voltage

x = Ha*((Hb*x1) + (B*(u(n,:)'+u(n-1,:)')) + (C*(I+I1)));% state update

y(n) = L*x + M*u(n,:)' + N*I;                           % find output

x1 = x;                                                 % step state
v1 = v;                                                 % step state
I1 = I;                                                 % step state

end

%%%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotting
    
    S = fft(y);                             % FFT of output
    S = S(1:floor(NT/2)) ;                  % trimmed above nyquist
    
    fvec = [0:floor(NT/2)-1]'*SR/NT ;       % frequency vector for plots
    
    figure(1)
    plot(tvec,u(:,1),tvec,y)
    legend('input','output','location','northeast')
    title({'Input and output waveforms from Range Master'})
    xlabel('Time (s)')
    ylabel('Magnitude')
    if input == "sin"
        xlim([0,0.05])
    end

    figure(2)
    loglog(fvec,abs(S))
    xlim([0,SR/2])
    title({'Frequency responce of output from Range Master'})
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    hold on
    
    figure(3)
    plot(tvec,u(:,2))
    title('supply voltage over time')
    xlabel('time (s)')
    ylabel('supply voltage (V)')
end

if sound
    
    soundsc(y,SR)                           % sound ouput

end

%   save sound
if SaveAudio
    y = y/max(y);
    filename = 'RangeMaster_Output_S1832740_Mason.wav';
    audiowrite(filename,y,SR);
    
end








