close all; clear all; clc
%Read POLARIS 
[~,p1_KM]=read_polaris('p110219_m125706_KM5KM6_1c_fvv0');
[~,p0_KM]=read_polaris('p110219_m125706_KM5KM6_0c_fvv0');
nav1_KM=read_nav('p110219_m125706_KM5KM6_1c_fvv0_nav');
nav0_KM=read_nav('p110219_m125706_KM5KM6_0c_fvv0_nav');

%***********************************************************************
%Set radar variables
c = physconst('LightSpeed');            %Speed of light             [m/s]
Tp = p0_KM.T;                           %Transmit pulse length      [s]
K = p0_KM.B/p0_KM.T;                    %Chirp rate                 [s^2]
fc = p0_KM.Fc;                          %Center frequency           [Hz]  
fs = p0_KM.Fs;                          %Sampling frequency         [Hz]
lambda = c/fc;                          %Wavelength                 [m]
PRF = p0_KM.PRF;                        %Pulse repetition frequency [Hz]
v = 75;                                 %Radar velocity             [m/s]

%Set simulation parameters
theta_s = (0*pi)/180;                   %Squint angle             
fdc = 2*v*sin(theta_s*pi/180)/lambda;   %Doppler centroid
Naz = 2048;                             %Raw data azimuth dim. size
Nra = 280;                              %Raw data range dim. size
BD = 50;                                %Doppler bandwidth of point targets     
La = 0.886*2*v*cos(theta_s)/BD;         %Real aperture length
beta_az = 0.886*lambda/La;              %Azimuth 3 dB beam width
%***********************************************************************

%Target locations (x0=azimuth,y0=range)
x0 = [0 -40 40 0 0];
y0 = [1200 1200 1200 1000 1400];
numTargs = length(x0);

%Populate target location matrix
pos = [x0;y0].';

%Determine the beam center crossing distances for the targets
for i = 1:numTargs
    xc(i) = x0(i)-y0(i)*tan(theta_s);   
end

%Range to the first point target, which is in the middle of the scene. 
Rc = y0(1);
%Closest approach range
R0 = Rc*cos(theta_s); 
%Synthetic aperture length
Ls = R0*beta_az;

%Imaging grid
gridxy = ones(Naz,Nra);
%Along-track matrix
x = repmat(v*((-Naz/2: Naz/2-1)/PRF),[Nra,1]).';
%Fast time matrix
y = repmat( 2*R0/c +(-Nra/2:Nra/2-1)/fs,[Naz,1]);
s0 = zeros(Naz,Nra);

%Start point-target simulation
for i = 1:numTargs   
    %Calculate range
    R = sqrt((y0(i).*gridxy).^2+(x-x0(i).*gridxy).^2);
    %Azimuth beam-pattern from small angle approx.
    wa= (sinc(0.886/beta_az .* atan( (x-xc(i))/y0(i) ) ) ).^2; 
    %Chirp envelope
    wr = (abs(y-2.*R/c)) <= ((Tp/2));
    %i'th received signal
    s_i = wr.*wa.*exp(-(1j*4*pi*fc).*R./c).*exp((1j*pi*K).*(y-2.*R./c).^2);
    %Accumulation of signal energy to the received data matrix
    s0 = s0 + s_i;
end

imagesc((abs(s0)))
figure
imagesc((real(s0)))
xlim([111 172]);ylim([708 1293])

%Time vector
t = -Tp/2:1/fs:Tp/2-1/fs;
s_Tx = exp(1i*pi*K*(t).^2);
%Zero padding reference signal
pad_length = Nra-length(s_Tx);
s_Tx_pad=[s_Tx zeros(1,pad_length)];

%Frequency vector for window function 
f_ra=[0:Nra-1]*fs/Nra-floor(([0:Nra-1]*fs/Nra)/(fs/2))*fs;
%Range weighting 
win_ra = 0.08+0.92*cos(pi*f_ra/p0_KM.B).^2;

%Perform range compression
for i=1:Naz 
y_ra = ifft((conj(fft(s_Tx_pad).*win_ra)).*fft(s0(i,:)));
y_RC(:,i) = y_ra(1:pad_length+1);
end

azaxis = (0:Naz-1)./PRF;
raaxis = (0:Nra-1)./fs;
figure
imagesc(azaxis,raaxis,(abs(y_RC)));xlabel('Along-track time');ylabel('Delay');

%***Perform backprojection focusing***

%Set range compressed data as input
y_in = y_RC;

%Set output image pixels to be focused
BP_ra = 1:size(y_RC,1)-3;
BP_az = 400:Naz-400; 

%Range image space
Del_R = y(1,:);
Del_RX =  Del_R(1);
rho_ra = c/2/fs;

%BP focusing
for i = BP_ra 
      R = c*Del_R(i)/2;  
      L = lambda*R/(2*(2.5)); 
      Np = round(L/(v/PRF));
      Np = Np+mod(Np,2);
      num_pulses = length(-Np/2:Np/2);
      win = blackman(num_pulses);         
    for j=BP_az
         x0 = (-Np/2:Np/2).*v/PRF;
         idx_az = j-Np/2:j+Np/2; 
         acc = 0;   
         for k = 1:num_pulses           
              trt = 2/c*(sqrt(x0(k)^2+R^2));
              Rcell = round((trt-Del_RX)*c/2/rho_ra)+1;          
              signal = y_in(Rcell,idx_az(k));                                  
              acc = acc + win(k).*signal.*exp(1i*2*pi*fc*trt);    
          %    acc = acc + signal.*exp(1i*2*pi*fc*trt);   
         end         
        y_foc(i,j) = acc;     
    end
i
end


