
clear all;
close all;

for mm = 4  % for the mmth dataset
    format long e;

    %****** prior values (c) and parameter ranges (cmin, cmax) ******************
          %    f1   k1     k2    N/C1      kn      kd      f_N2O_nit   f_N2O_dni         
       amin = [1e-6  1e-6   1e-6   0.001   1e-2    1e-6      0.001       5e-3]; 
       amax = [0.1   5e-2   8e-4   0.95    8e-1    9e-2      0.5         0.6];


    %******** input data *********
    sheet = strcat('sheet',num2str(mm));
    display(sheet);
    data_1 = xlsread('data_template_oneT.xlsx',sheet,'B3:AA10');
    data_1 = xlsread('data_template_oneT.xlsx',sheet,'B3:AA10');
    data_2 = xlsread('data_template_oneT_initials.xlsx','sheet1','A2:C1000');
    ctotal = data_2(mm,1);
    ntotal = data_2(mm,2);
    N_to_C_ratio = ntotal/ctotal;

    day_15_NH4 = data_1(1:find(isnan(data_1(:,1)))-1,1);
    NH4_15 = data_1(1:find(isnan(data_1(:,2)))-1,2);
    day_15_NO3 = data_1(1:find(isnan(data_1(:,3)))-1,3);
    NO3_15 = data_1(1:find(isnan(data_1(:,4)))-1,4);
    day_15_N2O = data_1(:,5);
    N2O_15 = data_1(:,6);
    day_15_CO2 = data_1(:,7);
    CO2_15 = data_1(:,8);

    %simulation starts
    adif  = (amax-amin)';
    a_op=amin'+rand*adif;
    a_new=zeros(8,1);      % create a new array with zeros  
    J_last = 300000 ;    % starting value for last 

    record_index=1;

    %Simulation starts
    upgraded=0;
    nsim    = 100000;      % number of simulations
    Nt = max([day_15_NH4(end),day_15_NO3(end),...
        day_15_N2O(end),day_15_CO2(end)]);
    for simu=1:nsim
        counter=simu
        upgraded % shows the upgraded number
        a_new=Generate_oneT(a_op,amin,amax);  %generate a new point

        f1 = a_new(1); 
        k1 = a_new(2);
        k2 = a_new(3);
        N_to_C1 = a_new(4);
        f2 = 1-f1;
        ntotal_1 = N_to_C1*f1*ctotal;
        ntotal_2 = ntotal - ntotal_1;
        N_to_C2 = ntotal_2/(f2*ctotal);

        kn = a_new(5);
        kd = a_new(6);
        
        f_N2O_nit = a_new(7);
        f_N2O_dni = a_new(8);

        Mod_NH4_15_last = NH4_15(1);
        Mod_NO3_15_last = NO3_15(1);
        
        for i = 1:Nt
           Mod_CO2_15_1(i) = (k1)*ctotal*f1*exp(-k1*i);
           Mod_CO2_15_2(i) = (k2)*ctotal*f2*exp(-k2*i);
           
           Mod_Nmin_15_1(i) = (k1*N_to_C1)*ntotal_1*exp(-k1*N_to_C1*i);
           Mod_Nmin_15_2(i) = (k2*N_to_C2)*ntotal_2*exp(-k2*N_to_C2*i);
           Mod_Nmin_15(i) =Mod_Nmin_15_1(i) + Mod_Nmin_15_2(i);
           
           Mod_Nnit_15(i) = Mod_NH4_15_last*kn;
           Mod_Ndni_15(i) = Mod_NO3_15_last*kd;
          
           Mod_NH4_15_present = Mod_NH4_15_last + Mod_Nmin_15(i) - Mod_Nnit_15(i);
           Mod_NO3_15_present = Mod_NO3_15_last + Mod_Nnit_15(i) - Mod_Ndni_15(i);
          
           Mod_NH4_15(i) = Mod_NH4_15_present;
           Mod_NO3_15(i) = Mod_NO3_15_present;
          
           Mod_NH4_15_last = Mod_NH4_15_present;
           Mod_NO3_15_last = Mod_NO3_15_present;
           
        end
        Mod_CO2_15 = Mod_CO2_15_1+Mod_CO2_15_2;
        Mod_N2O_15 = f_N2O_nit*Mod_Nnit_15 + f_N2O_dni*Mod_Ndni_15;
        
        for j = 1:length(day_15_NH4)-1
            Mod_NH4_15_n(j) = Mod_NH4_15(day_15_NH4(j+1));
        end;
        
        for j = 1:length(day_15_NO3)-1
            Mod_NO3_15_n(j) = Mod_NO3_15(day_15_NO3(j+1));
        end;
        
        for j = 1:length(day_15_N2O)-1
            Mod_N2O_15_n(j) = Mod_N2O_15(day_15_N2O(j+1));
        end;
        
        for j = 1:length(day_15_CO2)-1
            Mod_CO2_15_n(j) = Mod_CO2_15(day_15_CO2(j+1));
        end;
        

        J(1)  =  (norm(Mod_NH4_15_n' - NH4_15(2:end)))^2;
        J(2)  =  (norm(Mod_NO3_15_n' - NO3_15(2:end)))^2; 
        J(3)  =  (norm(Mod_N2O_15_n' - N2O_15(2:end)))^2; 
        J(4)  =  (norm(Mod_CO2_15_n' - CO2_15(2:end)))^2; 
        
        DJ(1) = 2*var(NH4_15(2:end));
        DJ(2) = 2*var(NO3_15(2:end));
        DJ(3) = 2*var(N2O_15(2:end));
        DJ(4) = 2*var(CO2_15(2:end));
        
        J_new= sum(J./DJ);

        % M-H algorithm
        delta_J = J_new-J_last;

        if min(1, exp(-delta_J)) >rand
            a_op=a_new;
            J_last=J_new;
            upgraded = upgraded+1;
            a_upgraded(:,upgraded)=a_op;  
            Mod_NH4_15_record(:,upgraded) = Mod_NH4_15_n';
            Mod_NO3_15_record(:,upgraded) = Mod_NO3_15_n';
            Mod_N2O_15_record(:,upgraded) = Mod_N2O_15_n';
            Mod_CO2_15_record(:,upgraded) = Mod_CO2_15_n';
            
            J_record(:,upgraded) = J_last;
        end

    end

    %%%%%%%%%%%%%%%
    %%% figures %%%
    %%%%%%%%%%%%%%%

    F1=figure(1);
    figure(1); plot(J_record); % 
    saveas(F1,'Fig 1','tif');
    
    % parameter ranges
    F3=figure(3);
    HS = 1500;
    HH = 100;
    halflen = ceil(upgraded/2);
    figure(3);subplot(2,4,1);hist(a_upgraded(1,halflen:end),HH);title('f1');xlabel('Parameter range');ylabel('frequency');axis([amin(1) amax(1) 1 HS]);
    figure(3);subplot(2,4,2);hist(a_upgraded(2,halflen:end),HH);title('k1');xlabel('Parameter range');ylabel('frequency');axis([amin(2) amax(2) 1 HS]);
    figure(3);subplot(2,4,3);hist(a_upgraded(3,halflen:end),HH);title('k2');xlabel('Parameter range');ylabel('frequency');axis([amin(3) amax(3) 1 HS]);
    figure(3);subplot(2,4,4);hist(a_upgraded(4,halflen:end),HH);title('N/C1');xlabel('Parameter range');ylabel('frequency');axis([amin(4) amax(4) 1 HS]);
    figure(3);subplot(2,4,5);hist(a_upgraded(5,halflen:end),HH);title('kn');xlabel('Parameter range');ylabel('frequency');axis([amin(5) amax(5) 1 HS]);
    figure(3);subplot(2,4,6);hist(a_upgraded(6,halflen:end),HH);title('kd');xlabel('Parameter range');ylabel('frequency');axis([amin(6) amax(6) 1 HS]);
    figure(3);subplot(2,4,7);hist(a_upgraded(7,halflen:end),HH);title('f-N2O-nit');xlabel('Parameter range');ylabel('frequency');axis([amin(7) amax(7) 1 HS]);
    figure(3);subplot(2,4,8);hist(a_upgraded(8,halflen:end),HH);title('f-N2O-dni');xlabel('Parameter range');ylabel('frequency');axis([amin(8) amax(8) 1 HS]);
    saveas(F3,'Fig 3','tif');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% write the parameters and new data into excel files %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %write the parameters into excel files
    colnames = {'f1','k1', 'k2',...
        'N/C1', 'kn', 'kd', 'f-N2O-nit', 'f-N2O-dni'};
    xlswrite('parameters_oneT',colnames, sheet, 'A1');
    xlswrite('parameters_oneT',a_upgraded(:,halflen:end)', sheet, 'A2');

    % write resp_mod to excel
    dat_1(:,1)=day_15_NH4(2:end); 
    dat_1(:,2)= mean(Mod_NH4_15_record')';
    dat_2(:,1)=day_15_NO3(2:end); 
    dat_2(:,2)= mean(Mod_NO3_15_record')';
    dat_3(:,1)=day_15_N2O(2:end); 
    dat_3(:,2)= mean(Mod_N2O_15_record')';
    dat_4(:,1)=day_15_CO2(2:end); 
    dat_4(:,2)= mean(Mod_CO2_15_record')';
    colnames2 = {'15','day','NH4', 'day','NO3', 'day','N2O', 'day','CO2'};
    xlswrite('Modeled_output_oneT',colnames2, sheet, 'A1');
    xlswrite('Modeled_output_oneT',dat_1, sheet, 'B2');
    xlswrite('Modeled_output_oneT',dat_2, sheet, 'D2');
    xlswrite('Modeled_output_oneT',dat_3, sheet, 'F2');
    xlswrite('Modeled_output_oneT',dat_4, sheet, 'H2');
    
    F2 = figure(2);
    OBS = [NH4_15(2:end);NO3_15(2:end);N2O_15(2:end);CO2_15(2:end)];
    MOD = [dat_1(:,2);dat_2(:,2);dat_3(:,2);dat_4(:,2)];
    fitlm=fitlm(OBS,MOD);
    Rsquire = fitlm.Rsquared.Ordinary;
    Rs = ['R^2 = ' num2str(Rsquire)];
    slope = OBS\MOD;
    MOD_CAL = slope*OBS;
    figure(2);subplot(3,2,1);plot(NH4_15(2:end),dat_1(:,2),'bo'); hold on;
    figure(2);subplot(3,2,2);plot(NO3_15(2:end),dat_2(:,2),'bo'); hold on;
    figure(2);subplot(3,2,3);plot(N2O_15(2:end),dat_3(:,2),'bo'); hold on;
    figure(2);subplot(3,2,4);plot(CO2_15(2:end),dat_4(:,2),'bo'); hold on;
    figure(2);subplot(3,2,5);plot(OBS,MOD,'bo'); hold on;
    figure(2);subplot(3,2,5);plot(OBS,MOD_CAL,'b-'); text(min(OBS),max(MOD),Rs);hold on;
    figure(2);subplot(3,2,6);plot(OBS,'bo'); hold on;
	figure(2);subplot(3,2,6);plot(MOD,'r^'); 
    hold off;
    saveas(F2,'Fig 2','tif');
    

end;