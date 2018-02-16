%% Reset Everything

close all
clear all

%% Time / Hour Angle
tm = 0:1:1439;  %Increment minutes in one day
th = tm/60;     %Convert so that it increments fraction of hrs in on day.
h = 12 - th;    %Convert to difference between solar noob
H = 15*h;       %Convert to Hour Angle 

%% Declimantion Angle
%  23.45sin[(360/365)(n-81)]

% Summer Day: July 1st
n_sum = 181; %days
n_win = 372; %days
sum_declimation = 23.45.*sind((360/365).*(n_sum - 81));
win_declimation = 23.45.*sind((360/365).*(n_win - 81));
        
%% Altitude Angle 
% Altitude angle at solar noon, Bn - angle between the sun and the local
% horizon:
    % Bn = 90 - Lat - Dec.
    
    % Location:
        % Latitude: 35.484 deg.
        % Longitude: 120.672 deg. 
        lat = 35.484;
        
    % Sin(beta) = cos(Lat)cos(Dec)cos(H) + sin(Lat)sin(Dec)
    sum_coslat = cosd(lat);
    sum_cosDec = cosd(sum_declimation);
    sum_cosH = cosd(H);
    sum_sinlat = sind(lat);
    sum_sinDec = sind(sum_declimation);
    sum_sinbeta = ((sum_coslat.*sum_cosDec.*sum_cosH) +  sum_sinlat.*sum_sinDec);
    sum_altitude_angle = asind(sum_sinbeta);
    
    win_coslat = cosd(lat);
    win_cosDec = cosd(win_declimation);
    win_cosH = cosd(H);
    win_sinlat = sind(lat);
    win_sinDec = sind(win_declimation);
    win_sinbeta = (win_coslat.*win_cosDec.*win_cosH) +  win_sinlat.*win_sinDec;
    win_altitude_angle = asind(win_sinbeta);
    
    
        
    %sum_altitude = 90 - lat - sum_declimation
    %win_altitude = 90 - lat - win_declimation
    
    %% Azimuth Angle 
% The azimuth angle, Az, of the sun or a celetial body (St) is depending on
% it's declination, Dec, and altitude, Alt, and on the Lat of the observer.

% azimuth angle: 
    % sin(phi) = [cos(Dec)sin(H)]/cos(altitude_angle)
    
sum_sinH = sind(H);
sum_cosBeta = cosd(sum_altitude_angle);

sum_sinphi = ((sum_cosDec.*sum_sinH)./sum_cosBeta);
sum_azimuth_angle = asind(sum_sinphi);

win_sinH = sind(H);
win_cosBeta = cosd(win_altitude_angle);

win_sinphi = (win_cosDec.*win_sinH)./win_cosBeta;
win_azimuth_angle = asind(win_sinphi);

 for(i=1:length(H))
     if (sum_cosH(i) >= (tand(sum_declimation)/tand(lat)))
     else if(sum_azimuth_angle(i)<0)
        sum_azimuth_angle(i) = -180 - sum_azimuth_angle(i);
     else 
     sum_azimuth_angle(i) = 180 - sum_azimuth_angle(i);
         end
     end
 end
 
 for(i=1:length(H))
     if(sum_altitude_angle(i) < 0)
         sum_altitude_angle(i) = 0;
     end
 end
 
  for(i=1:length(H))
     if (win_cosH(i) >= (tand(win_declimation)/tand(lat)))
     else if(win_azimuth_angle(i)<0)
        win_azimuth_angle(i) = -180 - win_azimuth_angle(i);
     else 
     win_azimuth_angle(i) = 180 - win_azimuth_angle(i);
         end
     end
  end
  
   for(i=1:length(H))
     if(win_altitude_angle(i) < 0)
         win_altitude_angle(i) = 0;
     end
   end
 
figure(1)
plot(sum_azimuth_angle , sum_altitude_angle)
hold on
plot(win_azimuth_angle , win_altitude_angle)


%% Irradiance 
% Irradiance of atmosphere: 1377 w/m^2

sum_AM = 1./sin(sum_altitude_angle);
%win_AM = 1./sin(win_altitude_angle);

sum_I = 1377.*(0.7.^(sum_AM).^(0.678));
%win_I = 1377.*(0.7.^(win_AM).^(0.678));
figure(2)
plot(tm, sum_I)



  
    