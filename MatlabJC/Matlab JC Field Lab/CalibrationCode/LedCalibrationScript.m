% LED calibration

% JC 9/16/2015 

%% params

Udt_power = .017e-9; % (Watts = J/sec) Udt meter measure
spotR = .0025 ; % (meters) spot size radius
CollectingArea = 0.5e-12; % MOUSE m^2

%% LED spectra
cd /Volumes/lab/Experiments/Calibration/LED' 'Spectra/ 
load LED_spectrum

wl = [400:1:800]*10^-9 ; % (meters) wavelength of measured led spectrum

%% Photoreceptor spectral sensitivty fit with Baylor nomogram
WL = wl*10^6 ; % put into um 
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;
lambdaMax = 491;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
PhotonSensitivity = 10.^LogPhotonSensitivity;

%% calculate photo isomerization rates

% scale LED spectrum according to Udt power measure
wl_energy = (3*10^8)*(6.63*10^-34)./wl ; % (Joules/photon) energy per photon wl

LED_energy = LED_spectrum.*wl_energy ; % scale energy by relative power of LED at each wl 

LED_power = sum(LED_energy) * (wl(2)-wl(1)) ; % integral of LED energy spectra

LED_energy_calibrated = LED_energy * Udt_power/LED_power ; % scale relative to Udt measure

LED_spectrum_calibrated = LED_energy_calibrated./wl_energy ; 

% impact of LED on photorecptor
Photons = PhotonSensitivity*LED_spectrum_calibrated' ;

PhotonCatch = CollectingArea*Photons/(pi*spotR^2) 




