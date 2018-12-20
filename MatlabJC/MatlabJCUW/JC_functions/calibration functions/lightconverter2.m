function Light = lightconverter2()


% enter and set the values from the light meter
meter = input('B1, B2, B3, B4, G1, G2, G3, R1, R2, R3 meter values * e-3 (enter as vector, NaN if no reading)\n')

% enter receptor type
CellType = input('receptor type?','s')

% isomerizations of interest
isomerization(1) = 5000 ;
isomerization(2) = 2.25 ;
isomerization(3) = 1 ;

% collecting area and divide by the attenuation of the filter paper
if strncmpi(CellType,'rod',3)==1 ; % if it is a rod
    collectingArea = .5 ; % um^2/rod
else
    collectingArea = .37 ;
end

filterPatten_red = 1.34 ; % these numbers are from Felice D.
filterPatten_green = 1.85 ;
filterPatten_blue = 2.85 ;


% set up output matrix
Light{1} = NaN(30,7) ;
Light{1}(:,1) = [1 2 3 4 1 2 3 1 2 3 1 2 3 4 1 2 3 1 2 3 1 2 3 4 1 2 3 1 2 3] ; % column 1 are the led swithes
Light{1}(:,2) = [meter';meter';meter'] ; % column 2 are the meter readings
Light{1}(:,3) = [0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 4 4 4 4 4 4 4 4 4 4] ; % column 3 is the NDF

Light{2} = {'BLUE1 NDF 0','BLUE2 NDF 0', 'BLUE3 NDF 0', 'BLUE4 NDF 0', 'brightgreen1 NDF 0', 'brightgreen2 NDF 0', 'brightgreen3 NDF 0','red1 NDF 0', 'red2 NDF 0', 'red3 NDF 0',...
    'BLUE1 NDf 2','BLUE2 NDF 2', 'BLUE3 NDF 2', 'BLUE4 NDF 2', 'brightgreen1 NDF 2', 'brightgreen2 NDF 2', 'brightgreen3 NDF 2','red1 NDF 2', 'red2 NDF 2', 'red3 NDF 2', ...
    'BLUE1 NDF 4','BLUE2 NDF 4', 'BLUE3 NDF 4', 'BLUE4 NDF 4', 'brightgreen1 NDF 4', 'brightgreen2 NDF 4', 'brightgreen3 NDF 4','red1 NDF 4', 'red2 NDF 4', 'red3 NDF 4'} ;

Light{3} = CellType ; 

% run BackgroundIntensity function for each meter reading

for NDF = 0:2:2 ; % for each NDF
for a = 1:4 % for each possible blue LED meter reading...
    if isnan(meter(a)) == 0 ; % if there is a meter reading given (NaN was not used)
        Light{1}(a+NDF*5,4) = (BackgroundIntensity_JC('Blue',NDF,meter(a),CellType)  * collectingArea)/filterPatten_blue ; % LED =BLUE, 
    end
end
for aa = 5:7 % for each possible green meter reading
   if isnan(meter(aa)) == 0 ; % if there is a meter reading given (NaN was not used)
        Light{1}(aa+NDF*5,4) = (BackgroundIntensity_JC('Green',NDF,meter(aa),CellType) * collectingArea)/filterPatten_green ; % LED =brightgreen,
    end
end 
for aaa = 8:10 % for each possible red meter reading
   if isnan(meter(aaa)) == 0 ; % if there is a meter reading given (NaN was not used)
        Light{1}(aaa+NDF*5,4) = (BackgroundIntensity_JC('Red',NDF,meter(aaa),CellType) * collectingArea)/filterPatten_red ; % LED =brightgreen,
    end
end     
end % end NDF loop


Light{1}(21:30,4) = Light{1}(11:20,4)./100 ; % make photon output of LED w/ NDF 4 100 times smaller than NDF of 2  

Light{1}(:,5) = isomerization(1)./ Light{1}(:,4) ; % column 5 is the voltage need to get 5000 P*/sec/receptor

Light{1}(:,6) = isomerization(2)./ Light{1}(:,4) ; % column 6 is the voltage need to get 150P*/sec/receptor

Light{1}(:,7) = isomerization(3)./ Light{1}(:,4); % column 7 is the voltage need to get 2.5 P*/sec/receptor


% find the setting and voltages best to use for 5000, 150, and 2.5 

vnearone = abs(Light{1}(:,5:7) - 1) ; %subtract 1 from each voltage
[v, nearone] = min(vnearone) ;   % for each isom rate find the voltage that is closest to one

% display results
for aaaa = 1:length(nearone) ; % for each best setting 

disp(isomerization(aaaa)) % display the desired isomerization rate
disp(Light{2}(nearone(aaaa))) % display the best setting to use
disp(Light{1}(nearone(aaaa),aaaa+4)) % display the voltage needed
end



