function PhotonFlux = BackgroundIntensity()%This function will ask you for the calibrated value of the setting used, the%area of the spot, the NDF in the light path and the type of cell.  The function%will return a value for background light this corresponds to in equivalent 500 nm %photons.%%Created and pillaged from LEDCalibration by APS Sep. 03LEDType = input('Which LED are you using?. \n', 's');calib = input('What is the calibrated value for the setting you are using in units of 1e-3 uW?. \n');clear NDF;NDF = input('What is the neutral density you are using?. \n')SpotSize = input('Enter the spot size of the LED in microns squared. \n');CellType = input('What type of photoreceptor are you stimulating?. \n', 's');% load the spectra of the photoreceptorif strcmpi(CellType, 'rod');	load macaque-rod-spectra	PhotoreceptorSpectra = MacaqueRodQuantalSpectraLog;elseif strcmpi(CellType, 'lcone')	load long_cone_spec	PhotoreceptorSpectra = SalLconeSpectraLog;elseif strcmpi(CellType, 'mcone');	load medium_cone_spec	PhotoreceptorSpectra = SalMconeSpectraLog;elseif strcmpi(CellType, 'scone');	load short_cone_spec	PhotoreceptorSpectra = SalSconeSpectraLog;elseif strcmpi(CellType, 'uvcone');	load uv_cone	PhotoreceptorSpectra = SaluvconeSpectraLog;endStartWaveLen = 400;EndWaveLen = 720;StepWaveLen = 1;TempReceptorSpectra = InterpolateNByTwo(PhotoreceptorSpectra, StartWaveLen, EndWaveLen, StepWaveLen);ReceptorSpectra(:,1) = TempReceptorSpectra(:,1);ReceptorSpectra(:,2) = 10.^TempReceptorSpectra(:,2);% Get spectra for LEDload DARK.TXTdark = DARK;LEDQuantalSpectra = dark(:,1);if strcmpi(LEDType, 'red')	load red.txt    red = RED	LEDQuantalSpectra(:, 2) = red(:, 2) - dark(:, 2);elseif strcmpi(LEDType, 'green')	load GREEN.TXT    green = GREEN;	LEDQuantalSpectra(:, 2) = green(:, 2) - dark(:, 2);elseif strcmpi(LEDType, 'brightgreen')	load brightgreen	clear LEDQuantalSpectra	LEDQuantalSpectra(:, 1) = brightgreen(:, 1);	LEDQuantalSpectra(:, 2) = brightgreen(:, 2);elseif strcmpi(LEDType, 'orange')	load orange.txt	LEDQuantalSpectra(:, 2) = orange(:, 2) - dark(:, 2);elseif strcmpi(LEDType, 'blue')	load BLUE.TXT    blue = BLUE;	LEDQuantalSpectra(:, 2) = blue(:, 2) - dark(:, 2);elseif strcmpi(LEDType, 'watermaze')	load watermaze.txt	LEDQuantalSpectra(:, 2) = watermaze(:, 2) - dark(:, 2);end% resample led spectra so they are at same points as photoreceptor spectraResampleFact = 1;ResampledLEDSpectra(:, 1) = Decimate(LEDQuantalSpectra(:, 1), ResampleFact);ResampledLEDSpectra(:, 2) = Decimate(LEDQuantalSpectra(:, 2), ResampleFact);LEDSpectra = InterpolateNByTwo(ResampledLEDSpectra, StartWaveLen, EndWaveLen, StepWaveLen);% the spectral properties of each of the NDF filters for each LED% Load blue NDF spectra because I don't have the spectra on the other LEDs right now.load BlueNDF;NDFSpectra = ReceptorSpectra;if NDF == 0;	NDFSpectra(:, 2) = 1endif NDF == 2;	if strcmpi(LEDType, 'red')		fprintf(1, 'Warning - no red ndf spectrum, using blue\n');		NDFSpectra = BlueNDFSpectra;		NDFSpectra(:,2) = NDFSpectra(:,2) * 0.01;			end	if strcmpi(LEDType, 'green')		fprintf(1, 'Warning - no green ndf spectrum, using blue\n');		NDFSpectra = BlueNDFSpectra;		NDFSpectra(:,2) = NDFSpectra(:,2) * 0.01;			end	if strcmpi(LEDType, 'brightgreen')		fprintf(1, 'Warning - no green ndf spectrum, using blue\n');		NDFSpectra = BlueNDFSpectra;		NDFSpectra(:,2) = NDFSpectra(:,2) * 0.01;			end	if strcmpi(LEDType, 'blue')		NDFSpectra = BlueNDFSpectra;		NDFSpectra(:,2) = NDFSpectra(:,2) * 0.01;	endendTempNDFSpectra = InterpolateNByTwo(NDFSpectra, StartWaveLen, EndWaveLen, StepWaveLen);clear NDFSpectraNDFSpectra = TempNDFSpectra;% --------------------------------------------------------------------------------------------------% Computes the photon flux%--------------------------------------------------------------------------------------------------			% convert quantal LED spectra to energy spectra and apply udt factors, then convert back to quantal spectra, % now in absolute unitsScFact = 6.6e-34 * 3e8;LEDEnergySpectra = LEDSpectra;LEDEnergySpectra(:, 2) = LEDSpectra(:, 2) ./ LEDSpectra(:, 1);LEDPow = sum(LEDEnergySpectra(:, 2)) * (LEDEnergySpectra(2, 1) - LEDEnergySpectra(1, 1));LEDEnergySpectra(:, 2) = LEDEnergySpectra(:, 2) * calib * 1e-9 / LEDPow;%1e-9 here puts the answer into WattsLEDSpectra(:, 2) = LEDEnergySpectra(:, 2) .* LEDEnergySpectra(:, 1) * 1e-9 / (ScFact * SpotSize);%1e-9 here converts wavelength to metersfiguresubplot(3,1,1)plot(LEDSpectra(:,1), LEDSpectra(:,2), 'b')title('LEDType')subplot(3,1,2)plot(ReceptorSpectra(:,1), ReceptorSpectra(:,2), 'r')subplot(3,1,3)plot(NDFSpectra(:,1), NDFSpectra(:,2), 'g')pause% multiply LED spectra by photoreceptor spectraPhotonFlux = sum(LEDSpectra(:, 2) .* ReceptorSpectra(:, 2) .* NDFSpectra(:, 2))