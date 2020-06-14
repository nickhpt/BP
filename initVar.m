function initVar(BP_ra,BP_az,scene)
disp('Initializing data...');
if contains(scene,'KM')
   [data,p] = read_polaris('p110219_m125706_KM5KM6_0c_fvv0');
   [~,p1] = read_polaris('p110219_m125706_KM5KM6_1c_fvv0');
   nav = read_nav('p110219_m125706_KM5KM6_1c_fvv0_nav'); 
   %Replace corrupt navigation data samples with neighbours
   nav.ve(1381)=nav.ve(1380);nav.ve(1382)=nav.ve(1380);
   nav.vn(1381)=nav.vn(1380);nav.vn(1382)=nav.vn(1380);
   nav.dif(1381)=nav.dif(1380);nav.dif(1382)=nav.dif(1380);
   nav.alt(1381)=nav.alt(1380);nav.alt(1382)=nav.alt(1380);
   
elseif contains(scene,'EP')
   [data,p]=read_polaris('p110218_m100213_epica_0c_svv0');
   [~,p1]=read_polaris('p110218_m100213_epica_1c_svv0');
    nav = read_nav('p110218_m100213_epica_1c_svv0_nav');
    
else %else initialize Jutulustraumen Glecier
    [data,p]=read_polaris('p110219_m112000_KM4KM5_0c_fvv0');
    [~,p1]=read_polaris('p110219_m112000_KM4KM5_1c_fvv0');
    nav = read_nav('p110219_m112000_KM4KM5_1c_fvv0_nav');

end


vars.PRF = p.PRF;
vars.t = (0:1/vars.PRF:p.Naz/vars.PRF-1/vars.PRF);
vars.c = physconst('LightSpeed');
vars.lambda = vars.c/p.Fc;          
vars.fc = p.Fc;
vars.Del_RX = p.RxDelay;
vars.fs = p.Fs;
vars.Del_R = (vars.Del_RX+(0:p.Nra-1)*1/vars.fs);
vars.rho_ra = 1/vars.fs*(vars.c/2);
vars.h_r = nav.dif;
vars.Zspacing = p1.ZSpacing;
vars.v = sqrt(nav.ve.^2+nav.vn.^2);
vars.Ra_cells=(((p1.ZOff-(nav.alt-nav.dif))/abs(...
vars.Zspacing)).');
vars.d_off = (vars.Ra_cells).*abs(vars.Zspacing);
vars.y_in = data;
vars.Azspacing = p1.AzSpacing;
vars.BP_ra = BP_ra;
vars.BP_az = BP_az;
save('datastruct','-struct','vars')


end

