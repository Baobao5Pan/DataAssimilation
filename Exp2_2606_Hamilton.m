
clear all;
close all;

for mm = 3  % for the mmth dataset
    format long e;

    %****** prior values (c) and parameter ranges (cmin, cmax) ******************
          %    f1  k1     k2    q1@15  q2@15  q1@25  q2@25  N/C1   kn    Qn   kd    Qd  f_N2O_nit   f_N2O_dni         
     a    = [0.05  1e-1  1e-5    2      2      2      2    0.10    1e-3   2   1e-3  2    0.005      0.005];   
   amin =  [0     0     1e-6    1      3      1       1    0.1     1e-5   1    1e-5  2    0          1e-4];    
   amax =  [0.10  0.1   3e-3    8      4      2.5      5    0.8    8e-2   6    5e-2  10   0.05        0.5];



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
    J_last = 100000 ;    % starting value for last 

    record_index=1;

    %Simulation starts
    upgraded=0;
    nsim    = 100000;      % number of simulations
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
%%
%%%%%%%%%%%%%%%
%%% Model Validation %%%
%%%%%%%%%%%%%%%    



frac_val = 0.2;
a_best = mean(a_upgraded(:, ceil(upgraded/2):end), 2);


val_idx_15_NH4 = floor(length(day_15_NH4)*(1-frac_val))+1:length(day_15_NH4);
val_idx_25_NH4 = floor(length(day_25_NH4)*(1-frac_val))+1:length(day_25_NH4);
val_idx_35_NH4 = floor(length(day_35_NH4)*(1-frac_val))+1:length(day_35_NH4);

val_idx_15_NO3 = floor(length(day_15_NO3)*(1-frac_val))+1:length(day_15_NO3);
val_idx_25_NO3 = floor(length(day_25_NO3)*(1-frac_val))+1:length(day_25_NO3);
val_idx_35_NO3 = floor(length(day_35_NO3)*(1-frac_val))+1:length(day_35_NO3);

val_idx_15_N2O = floor(length(day_15_N2O)*(1-frac_val))+1:length(day_15_N2O);
val_idx_25_N2O = floor(length(day_25_N2O)*(1-frac_val))+1:length(day_25_N2O);
val_idx_35_N2O = floor(length(day_35_N2O)*(1-frac_val))+1:length(day_35_N2O);

val_idx_15_CO2 = floor(length(day_15_CO2)*(1-frac_val))+1:length(day_15_CO2);
val_idx_25_CO2 = floor(length(day_25_CO2)*(1-frac_val))+1:length(day_25_CO2);
val_idx_35_CO2 = floor(length(day_35_CO2)*(1-frac_val))+1:length(day_35_CO2);


pred_NH4 = [Mod_NH4_15(day_15_NH4(val_idx_15_NH4)); ...
            Mod_NH4_25(day_25_NH4(val_idx_25_NH4)); ...
            Mod_NH4_35(day_35_NH4(val_idx_35_NH4))];

pred_NO3 = [Mod_NO3_15(day_15_NO3(val_idx_15_NO3)); ...
            Mod_NO3_25(day_25_NO3(val_idx_25_NO3)); ...
            Mod_NO3_35(day_35_NO3(val_idx_35_NO3))];

pred_N2O = [Mod_N2O_15(day_15_N2O(val_idx_15_N2O)); ...
            Mod_N2O_25(day_25_N2O(val_idx_25_N2O)); ...
            Mod_N2O_35(day_35_N2O(val_idx_35_N2O))];

pred_CO2 = [Mod_CO2_15(day_15_CO2(val_idx_15_CO2)); ...
            Mod_CO2_25(day_25_CO2(val_idx_25_CO2)); ...
            Mod_CO2_35(day_35_CO2(val_idx_35_CO2))];


obs_NH4 = [NH4_15(val_idx_15_NH4); NH4_25(val_idx_25_NH4); NH4_35(val_idx_35_NH4)];
obs_NO3 = [NO3_15(val_idx_15_NO3); NO3_25(val_idx_25_NO3); NO3_35(val_idx_35_NO3)];
obs_N2O = [N2O_15(val_idx_15_N2O); N2O_25(val_idx_25_N2O); N2O_35(val_idx_35_N2O)];
obs_CO2 = [CO2_15(val_idx_15_CO2); CO2_25(val_idx_25_CO2); CO2_35(val_idx_35_CO2)];


PRED_all = [pred_NH4(:); pred_NO3(:); pred_N2O(:); pred_CO2(:)];
OBS_all  = [obs_NH4(:);  obs_NO3(:);  obs_N2O(:);  obs_CO2(:)];


RMSE_all = sqrt(mean((PRED_all - OBS_all).^2));
R2_all = 1 - sum((PRED_all - OBS_all).^2) / sum((OBS_all - mean(OBS_all)).^2);

fprintf('\n======== Model Validation Result (All Temperatures) ========\n');
fprintf('RMSE (validation): %.4f\n', RMSE_all);
fprintf('R²   (validation): %.4f\n', R2_all);


figure;
scatter(OBS_all, PRED_all, 40, 'filled');
hold on;
plot([min(OBS_all) max(OBS_all)], [min(OBS_all) max(OBS_all)], 'r--', 'LineWidth', 1.5);
xlabel('Observed');
ylabel('Predicted');
title(sprintf('Validation: R² = %.3f, RMSE = %.3f', R2_all, RMSE_all));
grid on;

end;
%%
%%%%%%%%%%%%%%%
%%% Model uncertainty %%%
%%%%%%%%%%%%%%% 


OBS_all = [NH4_15(2:end); NO3_15(2:end); N2O_15(2:end); CO2_15(2:end); ...
           NH4_25(2:end); NO3_25(2:end); N2O_25(2:end); CO2_25(2:end); ...
           NH4_35(2:end); NO3_35(2:end); N2O_35(2:end); CO2_35(2:end)];

n_obs = length(OBS_all);
n_sim = size(Mod_NH4_15_record, 2);

SIM_all = zeros(n_obs, n_sim);
for i = 1:n_sim
    SIM_all(:,i) = [Mod_NH4_15_record(:,i); Mod_NO3_15_record(:,i); ...
                    Mod_N2O_15_record(:,i); Mod_CO2_15_record(:,i); ...
                    Mod_NH4_25_record(:,i); Mod_NO3_25_record(:,i); ...
                    Mod_N2O_25_record(:,i); Mod_CO2_25_record(:,i); ...
                    Mod_NH4_35_record(:,i); Mod_NO3_35_record(:,i); ...
                    Mod_N2O_35_record(:,i); Mod_CO2_35_record(:,i)];
end



SIM_median = median(SIM_all, 2);
SIM_lower  = prctile(SIM_all, 2.5, 2);
SIM_upper  = prctile(SIM_all, 97.5, 2);


figure;
hold on;
x = (1:n_obs)';
fill([x; flipud(x)], [SIM_lower; flipud(SIM_upper)], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(x, SIM_median, 'b-', 'LineWidth', 2);
plot(x, OBS_all, 'ko');
xlabel('Sampling points');
ylabel('Concentration or flux (mg/kg)');
legend('95% CI', 'Simulated', 'Observed');
title('Model Uncertainty Bounds');
grid on;

%%
%% =============================
%  Sensitivity Analysis on N2O (All Temperatures)
%  Each parameter ±20%, based on a_best
% =============================


a_best = mean(a_upgraded(:, ceil(upgraded/2):end), 2);

n_param = length(a_best);
S_index_up = zeros(n_param,1);
S_index_down = zeros(n_param,1);


Mod_N2O_base = run_N2O_allT(a_best, ctotal, ntotal, ...
    NH4_15(1), NH4_25(1), NH4_35(1), ...
    NO3_15(1), NO3_25(1), NO3_35(1), ...
    max([day_15_N2O; day_25_N2O; day_35_N2O]));

R_base = mean([mean(Mod_N2O_base{1}), mean(Mod_N2O_base{2}), mean(Mod_N2O_base{3})]);


for i = 1:n_param
    a_up = a_best;
    a_up(i) = a_best(i) * 1.2;

    a_down = a_best;
    a_down(i) = a_best(i) * 0.8;

    Mod_N2O_up = run_N2O_allT(a_up, ctotal, ntotal, ...
        NH4_15(1), NH4_25(1), NH4_35(1), ...
        NO3_15(1), NO3_25(1), NO3_35(1), ...
        max([day_15_N2O; day_25_N2O; day_35_N2O]));

    Mod_N2O_down = run_N2O_allT(a_down, ctotal, ntotal, ...
        NH4_15(1), NH4_25(1), NH4_35(1), ...
        NO3_15(1), NO3_25(1), NO3_35(1), ...
        max([day_15_N2O; day_25_N2O; day_35_N2O]));

    R_up = mean([mean(Mod_N2O_up{1}), mean(Mod_N2O_up{2}), mean(Mod_N2O_up{3})]);
    R_down = mean([mean(Mod_N2O_down{1}), mean(Mod_N2O_down{2}), mean(Mod_N2O_down{3})]);

    S_index_up(i) = (R_up - R_base) / R_base;
    S_index_down(i) = (R_down - R_base) / R_base;
end

figure;
bar([S_index_down, S_index_up]*100);
set(gca, 'XTickLabel', {'f1','k1','k2','Q1@15','Q2@15','Q1@25','Q2@25',...
                        'N/C1','kn','Qn','kd','Qd','fN2O_nit','fN2O_dni'});
legend('-20%', '+20%');
ylabel('Sensitivity Index (%)');
title('Parameter Sensitivity on Mean N₂O Flux (All Temperatures)');
xtickangle(45);

function Mod_N2O_all = run_N2O_allT(a, ctotal, ntotal, ...
    NH4_15_init, NH4_25_init, NH4_35_init, ...
    NO3_15_init, NO3_25_init, NO3_35_init, Nt)

    
    f1 = a(1); k1 = a(2); k2 = a(3);
    Q1_15 = a(4); Q2_15 = a(5);
    Q1_25 = a(6); Q2_25 = a(7);
    N_to_C1 = a(8);
    kn = a(9); Qn = a(10);
    kd = a(11); Qd = a(12);
    f_N2O_nit = a(13); f_N2O_dni = a(14);

    f2 = 1 - f1;
    ntotal_1 = N_to_C1*f1*ctotal;
    ntotal_2 = ntotal - ntotal_1;
    N_to_C2 = ntotal_2 / (f2*ctotal);

  
    Mod_N2O_15 = zeros(Nt,1);
    Mod_N2O_25 = zeros(Nt,1);
    Mod_N2O_35 = zeros(Nt,1);
    NH4_last_15 = NH4_15_init;
    NH4_last_25 = NH4_25_init;
    NH4_last_35 = NH4_35_init;
    NO3_last_15 = NO3_15_init;
    NO3_last_25 = NO3_25_init;
    NO3_last_35 = NO3_35_init;

    for t = 1:Nt
        % N mineralization
        Nmin_15 = k1/Q1_15*N_to_C1*ntotal_1*exp(-k1/Q1_15*N_to_C1*t) + ...
                  k2/Q2_15*N_to_C2*ntotal_2*exp(-k2/Q2_15*N_to_C2*t);
        Nmin_25 = k1*N_to_C1*ntotal_1*exp(-k1*N_to_C1*t) + ...
                  k2*N_to_C2*ntotal_2*exp(-k2*N_to_C2*t);
        Nmin_35 = k1*Q1_25*N_to_C1*ntotal_1*exp(-k1*Q1_25*N_to_C1*t) + ...
                  k2*Q2_25*N_to_C2*ntotal_2*exp(-k2*Q2_25*N_to_C2*t);

        % N transformation
        Nnit_15 = NH4_last_15 * kn/Qn;
        Nnit_25 = NH4_last_25 * kn;
        Nnit_35 = NH4_last_35 * kn*Qn;

        Ndni_15 = NO3_last_15 * kd/Qd;
        Ndni_25 = NO3_last_25 * kd;
        Ndni_35 = NO3_last_35 * kd*Qd;

        % Update NH4 and NO3
        NH4_present_15 = NH4_last_15 + Nmin_15 - Nnit_15;
        NO3_present_15 = NO3_last_15 + Nnit_15 - Ndni_15;
        NH4_present_25 = NH4_last_25 + Nmin_25 - Nnit_25;
        NO3_present_25 = NO3_last_25 + Nnit_25 - Ndni_25;
        NH4_present_35 = NH4_last_35 + Nmin_35 - Nnit_35;
        NO3_present_35 = NO3_last_35 + Nnit_35 - Ndni_35;

        % N2O simulation
        Mod_N2O_15(t) = f_N2O_nit*Nnit_15 + f_N2O_dni*Ndni_15;
        Mod_N2O_25(t) = f_N2O_nit*Nnit_25 + f_N2O_dni*Ndni_25;
        Mod_N2O_35(t) = f_N2O_nit*Nnit_35 + f_N2O_dni*Ndni_35;

        % Store for next iteration
        NH4_last_15 = NH4_present_15; NO3_last_15 = NO3_present_15;
        NH4_last_25 = NH4_present_25; NO3_last_25 = NO3_present_25;
        NH4_last_35 = NH4_present_35; NO3_last_35 = NO3_present_35;
    end

   
    Mod_N2O_all = {Mod_N2O_15, Mod_N2O_25, Mod_N2O_35};
end
