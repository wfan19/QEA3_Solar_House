function Angle = SunAngle(day_of_year, hr)
    %Sunrise/Sunset calculations
    %https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
    longitude = 42.3601 * pi / 180; %for Boston
    lat = -71.0589;
    timezone = -5;
    
    decl = -23.44 * cosd(360/365 * (day_of_year + 10));
    gamma = (2 * pi / 365) * (day_of_year - 1 + ((hr-12)/24));
    eqtime = 229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077 * sin(gamma) - 0.014615 * cos(2*gamma) - 0.040849*sin(2*gamma));
    time_offset = eqtime + 4*longitude - 60*timezone;
    tst = hr*60 + time_offset;
    ha = ((tst / 4) - 180); 
    solar_zenith = acosd((sind(lat)*sind(decl))+(cosd(lat)*cosd(decl)*cosd(ha)));
    solar_azimuth = 180 - acosd(-(sind(lat)*cosd(solar_zenith) - sind(decl))/(cosd(lat)*sind(solar_zenith)));
    Angle = [solar_zenith solar_azimuth];
end