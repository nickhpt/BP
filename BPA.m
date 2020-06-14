function [y_foc] = BPA(BP_ra,BP_az,scene)

initVar(BP_ra,BP_az,scene);
load('datastruct.mat');

M = 8;
hwait = waitbar(0,'Calculating...');
for i = BP_ra
      waitbar(i/BP_ra(end),hwait);
      R = c*Del_R(i)/2;  
      L = lambda*R/(2*(1.75));
      Np = round(L/(mean(v)/PRF));
      Np = Np+mod(Np,2);
      num_pulses = length(-Np/2:Np/2);
      win = blackman(num_pulses);
            
    for j=BP_az
         x0 = (-Np/2:Np/2).*v(j)/PRF;
         slowtimeBin = getSlowtimebins(j,mean(v),t);
         idx_az = slowtimeBin-Np/2:slowtimeBin+Np/2; 

        if i < Ra_cells(j); eps_r = 1;   else; eps_r = 3.15;  end     
   %    if i < Ra_cells(j); eps_r = 1;   else; eps_r = 1.69;  end  
         d = abs(Zspacing)*(i-1)-d_off(j);
         h = h_r(j);
         
         acc = 0;   
         for k = 1:num_pulses
              %Calculate round-trip delay
              r_air = h*sqrt(1+ ((x0(k))/( h+d/sqrt(eps_r) )  )^2  );    
              r_ice = d*sqrt(eps_r + ((x0(k))/(h+d/sqrt(eps_r)) )^2  );   
              trt = 2/c*(r_air+r_ice);   
              
              %Convert to imaging cells
              Rcell = round(fs*trt - fs*Del_RX + 1);
              Rcell2 = (fs*trt - fs*Del_RX + 1); 
              
              %Obtain input data (nearest neighbour)
              signal = y_in(Rcell,idx_az(k)); 
              
              %Obtain input data (sinc interpolator)
%               tidx = [floor(Rcell2)-(M-1):floor(Rcell2-1) floor(Rcell2)...
%                     ceil(Rcell2) ceil(Rcell2)+1:ceil(Rcell2)+(M-1)];
%               x = (y_in(tidx,idx_az(k)));
%               signal = sincInterp(x,tidx,Rcell2,M);  
              
              %Temporay phase and range history
%               phase1(k) = 2*pi*fc*(trt); 
%               Rcells(k) = Rcell; 
              
              %Accumulate weighted and phase corrected aperture
              acc = acc + win(k).*signal.*exp(1i*2*pi*fc*(trt)); 
         %    acc = acc + win(k).*signal; 
          %   acc = acc + signal.*exp(1i*2*pi*fc*(trt)); 
             
             
         end
         
        %Focus one pixel 
        y_foc(i,j) = acc;
                   
    end
    i
    %For each range save azimuth information for analysis
%     phases{i} = phase1;       %Save phase history 
%     azimuths{i} = x0;         %Save integration space 
%     vels(i) = v(j);           %Save velocities     
    
end

y_foc = y_foc(BP_ra,BP_az);
close(hwait)

end

