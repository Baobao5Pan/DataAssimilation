
clear all;
close all;

for mm = 1  % for the mmth dataset
    format long e;

    %****** prior values (c) and parameter ranges (cmin, cmax) ******************
          %    f1  k1     k2    q1@15  q2@15  q1@25  q2@25  N/C1   kn    Qn   kd    Qd  f_N2O_nit   f_N2O_dni         
     a    = [0.05  1e-1  1e-5    2      2      2      2    0.10    1e-3   2   1e-3  2    0.005      0.005];   
    amin =  [1e-6  1e-6  1e-5    1      2      1      1    0.1     1e-4   0.5  1e-5  1    0          1e-6];    
    amax =  [0.10  0.1   2e-3    3      6      5      3    0.8     2e-2   2    2e-2  6   1e-2        0.1];

    %******** input data *********
    sheet = strcat('sheet',num2str(mm));
    data_1 = xlsread('data_template_2.xlsx',sheet,'B3:AA10');
    data_2 = xlsread('data_template_2_initials.xlsx','sheet1','A2:C1000');
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

    day_25_NH4 = data_1(1:find(isnan(data_1(:,10)))-1,10);
    NH4_25 = data_1(1:find(isnan(data_1(:,11)))-1,11);
    day_25_NO3 = data_1(1:find(isnan(data_1(:,12)))-1,12);
    NO3_25 = data_1(1:find(isnan(data_1(:,13)))-1,13);
    day_25_N2O = data_1(:,14);
    N2O_25 = data_1(:,15);
    day_25_CO2 = data_1(:,16);
    CO2_25 = data_1(:,17);

    day_35_NH4 = data_1(1:find(isnan(data_1(:,19)))-1,19);
    NH4_35 = data_1(1:find(isnan(data_1(:,20)))-1,20);
    day_35_NO3 = data_1(1:find(isnan(data_1(:,21)))-1,21);
    NO3_35 = data_1(1:find(isnan(data_1(:,22)))-1,22);
    day_35_N2O = data_1(:,23);
    N2O_35 = data_1(:,24);
    day_35_CO2 = data_1(:,25);
    CO2_35 = data_1(:,26);

    %simulation starts
    adif  = (amax-amin)';
    a_op=amin'+rand*adif;
    a_new=zeros(14,1);      % create a new array with zeros  
    J_last = 300000 ;    % starting value for last 

    record_index=1;

    %Simulation starts
    upgraded=0;
    nsim    = 300000;      % number of simulations
    Nt = max([day_15_NH4(end),day_25_NH4(end),day_35_NH4(end),day_15_NO3(end),day_25_NO3(end),day_35_NO3(end),...
        day_15_N2O(end),day_25_N2O(end),day_35_N2O(end),day_15_CO2(end),day_25_CO2(end),day_35_CO2(end)]);
    for simu=1:nsim
        counter=simu
        upgraded % shows the upgraded number
        a_new=Generate_1(a_op,amin,amax);  %generate a new point

        f1 = a_new(1); 
        k1 = a_new(2);
        k2 = a_new(3);
        Q1_15 = a_new(4);   
        Q2_15 = a_new(5);
        Q1_25 = a_new(6);   
        Q2_25 = a_new(7);
        N_to_C1 = a_new(8);
        f2 = 1-f1;
        ntotal_1 = N_to_C1*f1*ctotal;
        ntotal_2 = ntotal - ntotal_1;
        N_to_C2 = ntotal_2/(f2*ctotal);

        kn = a_new(9);
        Qn = a_new(10);
        kd = a_new(11);
        Qd = a_new(12);

        f_N2O_nit = a_new(13);
        f_N2O_dni = a_new(14);

        Mod_NH4_15_last = NH4_15(1);
        Mod_NH4_25_last = NH4_25(1);
        Mod_NH4_35_last = NH4_35(1);
        Mod_NO3_15_last = NO3_15(1);
        Mod_NO3_25_last = NO3_25(1);
        Mod_NO3_35_last = NO3_35(1);

        for i = 1:Nt
           Mod_CO2_15_1(i) = (k1/Q1_15)*ctotal*f1*exp(-k1/Q1_15*i);
           Mod_CO2_15_2(i) = (k2/Q2_15)*ctotal*f2*exp(-k2/Q2_15*i);
           Mod_CO2_25_1(i) = (k1)*ctotal*f1*exp(-k1*i);
           Mod_CO2_25_2(i) = (k2)*ctotal*f2*exp(-k2*i);
           Mod_CO2_35_1(i) = (k1*Q1_25)*ctotal*f1*exp(-k1*Q1_25*i);
           Mod_CO2_35_2(i) = (k2*Q2_25)*ctotal*f2*exp(-k2*Q2_25*i);

           Mod_Nmin_15_1(i) = (k1/Q1_15*N_to_C1)*ntotal_1*exp(-k1/Q1_15*N_to_C1*i);
           Mod_Nmin_15_2(i) = (k2/Q2_15*N_to_C2)*ntotal_2*exp(-k2/Q2_15*N_to_C2*i);
           Mod_Nmin_15(i) =Mod_Nmin_15_1(i) + Mod_Nmin_15_2(i);
           Mod_Nmin_25_1(i) = (k1*N_to_C1)*ntotal_1*exp(-k1*N_to_C1*i);
           Mod_Nmin_25_2(i) = (k2*N_to_C2)*ntotal_2*exp(-k2*N_to_C2*i);
           Mod_Nmin_25(i) =Mod_Nmin_25_1(i) + Mod_Nmin_25_2(i);
           Mod_Nmin_35_1(i) = (k1*Q1_25*N_to_C1)*ntotal_1*exp(-k1*Q1_25*N_to_C1*i);
           Mod_Nmin_35_2(i) = (k2*Q2_25*N_to_C2)*ntotal_2*exp(-k2*Q2_25*N_to_C2*i);
           Mod_Nmin_35(i) =Mod_Nmin_35_1(i) + Mod_Nmin_35_2(i);

           Mod_Nnit_15(i) = Mod_NH4_15_last*kn/Qn;
           Mod_Nnit_25(i) = Mod_NH4_25_last*kn;
           Mod_Nnit_35(i) = Mod_NH4_35_last*kn*Qn;
           Mod_Ndni_15(i) = Mod_NO3_15_last*kd/Qd;
           Mod_Ndni_25(i) = Mod_NO3_25_last*kd;
           Mod_Ndni_35(i) = Mod_NO3_35_last*kd*Qd;

         
           Mod_NH4_15_present = Mod_NH4_15_last + Mod_Nmin_15(i) - Mod_Nnit_15(i);
           Mod_NH4_25_present = Mod_NH4_25_last + Mod_Nmin_25(i) - Mod_Nnit_15(i);
           Mod_NH4_35_present = Mod_NH4_35_last + Mod_Nmin_35(i) - Mod_Nnit_15(i);
           Mod_NO3_15_present = Mod_NO3_15_last + Mod_Nnit_15(i) - Mod_Ndni_15(i);
           Mod_NO3_25_present = Mod_NO3_25_last + Mod_Nnit_25(i) - Mod_Ndni_25(i);
           Mod_NO3_35_present = Mod_NO3_35_last + Mod_Nnit_35(i) - Mod_Ndni_35(i);

           Mod_NH4_15(i) = Mod_NH4_15_present;
           Mod_NH4_25(i) = Mod_NH4_25_present;
           Mod_NH4_35(i) = Mod_NH4_35_present;
           Mod_NO3_15(i) = Mod_NO3_15_present;
           Mod_NO3_25(i) = Mod_NO3_25_present;
           Mod_NO3_35(i) = Mod_NO3_35_present;

           Mod_NH4_15_last = Mod_NH4_15_present;
           Mod_NH4_25_last = Mod_NH4_25_present;
           Mod_NH4_35_last = Mod_NH4_35_present;
           Mod_NO3_15_last = Mod_NO3_15_present;
           Mod_NO3_25_last = Mod_NO3_25_present;
           Mod_NO3_35_last = Mod_NO3_35_present;

        end
        Mod_CO2_15 = Mod_CO2_15_1+Mod_CO2_15_2;
        Mod_CO2_25 = Mod_CO2_25_1+Mod_CO2_25_2;
        Mod_CO2_35 = Mod_CO2_35_1+Mod_CO2_35_2;
        Mod_N2O_15 = f_N2O_nit*Mod_Nnit_15 + f_N2O_dni*Mod_Ndni_15;
        Mod_N2O_25 = f_N2O_nit*Mod_Nnit_25 + f_N2O_dni*Mod_Ndni_25;
        Mod_N2O_35 = f_N2O_nit*Mod_Nnit_35 + f_N2O_dni*Mod_Ndni_35;

        for j = 1:length(day_15_NH4)-1
            Mod_NH4_15_n(j) = Mod_NH4_15(day_15_NH4(j+1));
        end;
        for j = 1:length(day_25_NH4)-1
            Mod_NH4_25_n(j) = Mod_NH4_25(day_25_NH4(j+1));
        end;
        for j = 1:length(day_35_NH4)-1
            Mod_NH4_35_n(j) = Mod_NH4_35(day_35_NH4(j+1));
        end;

        for j = 1:length(day_15_NO3)-1
            Mod_NO3_15_n(j) = Mod_NO3_15(day_15_NO3(j+1));
        end;
        for j = 1:length(day_25_NO3)-1
            Mod_NO3_25_n(j) = Mod_NO3_25(day_25_NO3(j+1));
        end;
        for j = 1:length(day_35_NO3)-1
            Mod_NO3_35_n(j) = Mod_NO3_35(day_35_NO3(j+1));
        end;

        for j = 1:length(day_15_N2O)-1
            Mod_N2O_15_n(j) = Mod_N2O_15(day_15_N2O(j+1));
        end;
        for j = 1:length(day_25_N2O)-1
            Mod_N2O_25_n(j) = Mod_N2O_25(day_25_N2O(j+1));
        end;
        for j = 1:length(day_35_N2O)-1
            Mod_N2O_35_n(j) = Mod_N2O_35(day_35_N2O(j+1));
        end;

        for j = 1:length(day_15_CO2)-1
            Mod_CO2_15_n(j) = Mod_CO2_15(day_15_CO2(j+1));
        end;
        for j = 1:length(day_25_CO2)-1
            Mod_CO2_25_n(j) = Mod_CO2_25(day_25_CO2(j+1));
        end;
        for j = 1:length(day_35_CO2)-1
            Mod_CO2_35_n(j) = Mod_CO2_35(day_35_CO2(j+1));
        end;


        J(1)  =  (norm(Mod_NH4_15_n' - NH4_15(2:end)))^2;
        J(2)  =  (norm(Mod_NO3_15_n' - NO3_15(2:end)))^2; 
        J(3)  =  (norm(Mod_N2O_15_n' - N2O_15(2:end)))^2; 
        J(4)  =  (norm(Mod_CO2_15_n' - CO2_15(2:end)))^2; 
        J(5)  =  (norm(Mod_NH4_25_n' - NH4_25(2:end)))^2;
        J(6)  =  (norm(Mod_NO3_25_n' - NO3_25(2:end)))^2; 
        J(7)  =  (norm(Mod_N2O_25_n' - N2O_25(2:end)))^2; 
        J(8)  =  (norm(Mod_CO2_25_n' - CO2_25(2:end)))^2; 
        J(9)  =  (norm(Mod_NH4_35_n' - NH4_35(2:end)))^2;
        J(10)  =  (norm(Mod_NO3_35_n' - NO3_35(2:end)))^2; 
        J(11)  =  (norm(Mod_N2O_35_n' - N2O_35(2:end)))^2; 
        J(12)  =  (norm(Mod_CO2_35_n' - CO2_35(2:end)))^2; 

        DJ(1) = 2*var(NH4_15(2:end));
        DJ(2) = 2*var(NO3_15(2:end));
        DJ(3) = 2*var(N2O_15(2:end));
        DJ(4) = 2*var(CO2_15(2:end));
        DJ(5) = 2*var(NH4_25(2:end));
        DJ(6) = 2*var(NO3_25(2:end));
        DJ(7) = 2*var(N2O_25(2:end));
        DJ(8) = 2*var(CO2_25(2:end));
        DJ(9) = 2*var(NH4_35(2:end));
        DJ(10) = 2*var(NO3_35(2:end));
        DJ(11) = 2*var(N2O_35(2:end));
        DJ(12) = 2*var(CO2_35(2:end));

        J_new= sum(J./DJ);

        % M-H algorithm
        delta_J = J_new-J_last;

        if min(1, exp(-delta_J)) >rand
            a_op=a_new;
            J_last=J_new;
            upgraded = upgraded+1;
            a_upgraded(:,upgraded)=a_op;  
            Mod_NH4_15_record(:,upgraded) = Mod_NH4_15_n';
            Mod_NH4_25_record(:,upgraded) = Mod_NH4_25_n';
            Mod_NH4_35_record(:,upgraded) = Mod_NH4_35_n';
            Mod_NO3_15_record(:,upgraded) = Mod_NO3_15_n';
            Mod_NO3_25_record(:,upgraded) = Mod_NO3_25_n';
            Mod_NO3_35_record(:,upgraded) = Mod_NO3_35_n';
            Mod_N2O_15_record(:,upgraded) = Mod_N2O_15_n';
            Mod_N2O_25_record(:,upgraded) = Mod_N2O_25_n';
            Mod_N2O_35_record(:,upgraded) = Mod_N2O_35_n';
            Mod_CO2_15_record(:,upgraded) = Mod_CO2_15_n';
            Mod_CO2_25_record(:,upgraded) = Mod_CO2_25_n';
            Mod_CO2_35_record(:,upgraded) = Mod_CO2_35_n';

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
    HS = 500;
    HH = 100;
    halflen = ceil(upgraded/2);
    figure(3);subplot(3,5,1);hist(a_upgraded(1,halflen:end),HH);title('f1');xlabel('Parameter range');ylabel('frequency');axis([amin(1) amax(1) 1 HS]);
    figure(3);subplot(3,5,2);hist(a_upgraded(2,halflen:end),HH);title('k1');xlabel('Parameter range');ylabel('frequency');axis([amin(2) amax(2) 1 HS]);
    figure(3);subplot(3,5,3);hist(a_upgraded(3,halflen:end),HH);title('k2');xlabel('Parameter range');ylabel('frequency');axis([amin(3) amax(3) 1 HS]);
    figure(3);subplot(3,5,4);hist(a_upgraded(4,halflen:end),HH);title('Q1@15');xlabel('Parameter range');ylabel('frequency');axis([amin(4) amax(4) 1 HS]);
    figure(3);subplot(3,5,5);hist(a_upgraded(5,halflen:end),HH);title('Q2@15');xlabel('Parameter range');ylabel('frequency');axis([amin(5) amax(5) 1 HS]);
    figure(3);subplot(3,5,6);hist(a_upgraded(6,halflen:end),HH);title('Q1@25');xlabel('Parameter range');ylabel('frequency');axis([amin(6) amax(6) 1 HS]);
    figure(3);subplot(3,5,7);hist(a_upgraded(7,halflen:end),HH);title('Q2@25');xlabel('Parameter range');ylabel('frequency');axis([amin(7) amax(7) 1 HS]);
    figure(3);subplot(3,5,8);hist(a_upgraded(8,halflen:end),HH);title('N/C1');xlabel('Parameter range');ylabel('frequency');axis([amin(8) amax(8) 1 HS]);
    figure(3);subplot(3,5,9);hist(a_upgraded(9,halflen:end),HH);title('kn');xlabel('Parameter range');ylabel('frequency');axis([amin(9) amax(9) 1 HS]);
    figure(3);subplot(3,5,10);hist(a_upgraded(10,halflen:end),HH);title('Qn');xlabel('Parameter range');ylabel('frequency');axis([amin(10) amax(10) 1 HS]);
    figure(3);subplot(3,5,11);hist(a_upgraded(11,halflen:end),HH);title('kd');xlabel('Parameter range');ylabel('frequency');axis([amin(11) amax(11) 1 HS]);
    figure(3);subplot(3,5,12);hist(a_upgraded(12,halflen:end),HH);title('Qd');xlabel('Parameter range');ylabel('frequency');axis([amin(12) amax(12) 1 HS]);
    figure(3);subplot(3,5,13);hist(a_upgraded(13,halflen:end),HH);title('f-N2O-nit');xlabel('Parameter range');ylabel('frequency');axis([amin(13) amax(13) 1 HS]);
    figure(3);subplot(3,5,14);hist(a_upgraded(14,halflen:end),HH);title('f-N2O-dni');xlabel('Parameter range');ylabel('frequency');axis([amin(14) amax(14) 1 HS]);
    saveas(F3,'Fig 3','tif');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% write the parameters and new data into excel files %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %write the parameters into excel files
    colnames = {'f1','k1', 'k2', 'Q1@15', 'Q2@15', 'Q1@25', 'Q2@25',...
        'N/C1', 'kn', 'Qn', 'kd', 'Qd', 'f-N2O-nit', 'f-N2O-dni'};
    xlswrite('parameters',colnames, sheet, 'A1');
    xlswrite('parameters',a_upgraded(:,halflen:end)', sheet, 'A2');

    % write resp_mod to excel
    dat_1(:,1)=day_15_NH4(2:end); 
    dat_1(:,2)= mean(Mod_NH4_15_record')';
    dat_2(:,1)=day_15_NO3(2:end); 
    dat_2(:,2)= mean(Mod_NO3_15_record')';
    dat_3(:,1)=day_15_N2O(2:end); 
    dat_3(:,2)= mean(Mod_N2O_15_record')';
    dat_4(:,1)=day_15_CO2(2:end); 
    dat_4(:,2)= mean(Mod_CO2_15_record')';
    dat_5(:,1)=day_25_NH4(2:end); 
    dat_5(:,2)= mean(Mod_NH4_25_record')';
    dat_6(:,1)=day_25_NO3(2:end); 
    dat_6(:,2)= mean(Mod_NO3_25_record')';
    dat_7(:,1)=day_25_N2O(2:end); 
    dat_7(:,2)= mean(Mod_N2O_25_record')';
    dat_8(:,1)=day_25_CO2(2:end); 
    dat_8(:,2)= mean(Mod_CO2_25_record')';
    dat_9(:,1)=day_35_NH4(2:end); 
    dat_9(:,2)= mean(Mod_NH4_35_record')';
    dat_10(:,1)=day_35_NO3(2:end); 
    dat_10(:,2)= mean(Mod_NO3_35_record')';
    dat_11(:,1)=day_35_N2O(2:end); 
    dat_11(:,2)= mean(Mod_N2O_35_record')';
    dat_12(:,1)=day_35_CO2(2:end); 
    dat_12(:,2)= mean(Mod_CO2_35_record')';
    colnames2 = {'15','day','NH4', 'day','NO3', 'day','N2O', 'day','CO2',...
        '25','day','NH4', 'day','NO3', 'day','N2O', 'day','CO2',...
        '35','day','NH4', 'day','NO3', 'day','N2O', 'day','CO2'};
    xlswrite('Modeled_output',colnames2, sheet, 'A1');
    xlswrite('Modeled_output',dat_1, sheet, 'B2');
    xlswrite('Modeled_output',dat_2, sheet, 'D2');
    xlswrite('Modeled_output',dat_3, sheet, 'F2');
    xlswrite('Modeled_output',dat_4, sheet, 'H2');
    xlswrite('Modeled_output',dat_5, sheet, 'K2');
    xlswrite('Modeled_output',dat_6, sheet, 'M2');
    xlswrite('Modeled_output',dat_7, sheet, 'O2');
    xlswrite('Modeled_output',dat_8, sheet, 'Q2');
    xlswrite('Modeled_output',dat_9, sheet, 'T2');
    xlswrite('Modeled_output',dat_10, sheet, 'V2');
    xlswrite('Modeled_output',dat_11, sheet, 'X2');
    xlswrite('Modeled_output',dat_12, sheet, 'Z2');
    
   F2 = figure(2);
    OBS = [NH4_15(2:end);NO3_15(2:end); N2O_15(2:end);CO2_15(2:end);NH4_25(2:end);NO3_25(2:end); N2O_25(2:end);CO2_25(2:end);NH4_35(2:end);NO3_35(2:end); N2O_35(2:end);CO2_35(2:end)];
    MOD = [dat_1(:,2);dat_2(:,2);dat_3(:,2);dat_4(:,2);dat_5(:,2);  dat_6(:,2); dat_7(:,2); dat_8(:,2); dat_9(:,2); dat_10(:,2);dat_11(:,2); dat_12(:,2)];
    fitlm=fitlm(OBS,MOD);
    Rsquire = fitlm.Rsquared.Ordinary;
    Rs = ['R^2 = ' num2str(Rsquire)];
    slope = OBS\MOD;
    MOD_CAL = slope*OBS;
    
    figure(2);subplot(3,2,1);plot(OBS,MOD,'bo'); hold on;
    figure(2);subplot(3,2,1);plot(OBS,MOD_CAL,'b-'); text(min(OBS),max(MOD),Rs);hold on;
    figure(2);subplot(3,2,2);plot(OBS,'bo'); hold on;
	figure(2);subplot(3,2,2);plot(MOD,'r^'); 
    hold off;
    saveas(F2,'Fig 2','tif');
end;