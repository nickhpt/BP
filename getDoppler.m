function [theta_max,fdmax,tvec,dPHdt] = getDoppler(phases,azimuths,vels,BP_ra)

%Computes instantaneous Doppler and maximum look angle

%Inputs (for each range)
%1.phases: Phase history
%2.azimuths: Along-track space
%3.vels: velocities
%4.BP_ra: Output range dimension vector

%Outputs (for each range)
%1.theta_max: Maximum look angles
%2.fdmax: Maximum Doppler frequencies
%3.tvec: Time vectors for plotting
%4.dPHdt: Doppler frequency

phases2 = phases(BP_ra);
azimuths2 = azimuths(BP_ra);
vels2 = vels(BP_ra);
lambda = 0.6892;
N = length(phases2);

for p = 1:N
    tvec{1,p} = cell2mat(azimuths2(1,p))./vels2(p);                                
    dPHdt{1,p} = (diff(cell2mat(phases2(1,p)))./diff(cell2mat(tvec(1,p))))/(2*pi); 
    ttemp = tvec{1,p};
    tdop{1,p} = (ttemp(2:end)+ttemp(1:(end-1)))/2;                                
    fdmax(p) = max(dPHdt{1,p});                                                    
    theta_max(p) = asin(fdmax(p)*lambda/(2*(vels2(p))))*180/pi;              
end

end

