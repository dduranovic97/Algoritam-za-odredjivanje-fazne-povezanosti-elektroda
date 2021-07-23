%------------------------------------------------------------------------%
%                  Daniel Duranovic, Natalija Ivos                       %
%           Automatizirana instrumentacija - Seminar (tema 3.)           %
%                         Akad.god. 2020./21.                            %
%------------------------------------------------------------------------%
%% Napomena: Radeno u MatLabu R2018a
%% Setup
clc
clear all
close all
%% Load data
% EEG - struktura koja sadrži 16 kanalni EEG i podatke o mjerenju
% EEG.times - vektor - vrijeme u ms
% EEG.srate - sampling rate (512 - fs)

load EEG_otvorene_oci.mat 
% load EEG_zatvorene_oci.mat % u?itati jedan od dva file-a

fs = EEG.srate;
N = fs; 
freqs = [0:511];
time = EEG.times/1000;

data = EEG.data; %16 EEG kanala
names = EEG.chaninfo.filecontent; %sadrži imena kanala

%% Uzima pozicije i imena elektroda i sprema u matricu (lakša usporedba koordinata i imena)
elektroda = {};
for i = 1:length(names(:,1))
    elektroda{i} = textscan(names(i,:), '%d %f %f %s');%[broj, kut, magn, ime]
end

%% Filtracija (PP) na 10 Hz
data_filt = [];
Fpass = [9.99,10.01];
for i = 1:length(data(:,1))
    x = bandpass(data(i,:), Fpass, fs);
    data_filt = [data_filt;x];
end

%% Usporedba (logika) - zanemaruje istu grupu 
%sprema u nxn matricu (n = broj elektroda)

MSC_10Hz = NaN(length(data(:,1))); %nxn NaN matrica
f_msc = 10; %frekv na kojoj uzimamo vr MSC-a
MaxMSC = [0,0,0]; %[elek1, elek2, vr]
MinMSC = [0,0,1];

for i = 1:length(elektroda)
    for j = 1:length(elektroda)
        
        if isnan(MSC_10Hz(i,j)) == 0 %da ne racuna 2 puta iste elektrode
            continue
        
            else if strncmpi(elektroda{i}{4},elektroda{j}{4}, 1) %ista zona (po imenu)
                continue
         
                    else 
                        %algoritam
                        dx = mscohere(data_filt(i,:),data_filt(j,:),hann(fs),[],[f_msc-1,f_msc+1],fs); %window, noverlap, [fmin,fmax], fs
              
                        MSC_10Hz(i,j) = dx(2); %10 Hz
                        MSC_10Hz(j,i) = dx(2);

                        fprintf('--Usporedba-- %s i %s; MSC: %f \r', char(elektroda{i}{4}),char(elektroda{j}{4}),dx(2))

                        %Min. i max vr.
                        if dx(2) > MaxMSC(3) 
                            MaxMSC(3) = dx(2);
                            MaxMSC(1) = i;
                            MaxMSC(2) = j;

                            else if dx(2) < MinMSC(3) 
                                    MinMSC(3) = dx(2);
                                    MinMSC(1) = i;
                                    MinMSC(2) = j;
                                end
                        end
            end %end if zona
        end %end da ne racuna 2 puta istu
    end %end j-petlje
end %end i-petlje

%%  Usporedjuje sve elektrode (treba uncomentat pa zakomentirati ovaj iznad)
% for i = 1:length(elektroda)
%     for j = 1:length(elektroda)
%         
%         %da ne racuna 2 puta iste elektrode
%         if isnan(MSC_10Hz(i,j)) == 0 
%             continue
%         end
%          
%             %algoritam
%             dx = mscohere(data_filt(i,:),data_filt(j,:),hann(fs),[],[f_msc-1,f_msc+1],fs); %window, noverlap, [fmin,fmax], fs
%             MSC_10Hz(i,j) = dx(2);
%             MSC_10Hz(j,i) = dx(2);
%             
%             fprintf('--Usporedba-- %d i %d; MSC: %f \r', elektroda{i}{1},elektroda{j}{1},dx(2))
%             
%             %Min. i max vr.
%             if dx(2) > MaxMSC(3) 
%                 MaxMSC(3) = dx(2);
%                 MaxMSC(1) = i;
%                 MaxMSC(2) = j;
%                 
%             else if dx(2) < MinMSC(3) 
%                     MinMSC(3) = dx(2);
%                     MinMSC(1) = i;
%                     MinMSC(2) = j;
%                 end
%             end
%         end
% end

%% Min i Max MSC [elektroda1, elektroda2, vrijednost]
fprintf('Elektroda %d i elektroda %d imaju najbolju povezanost (MSC :%f)\r', MaxMSC(1),MaxMSC(2),MaxMSC(3)); 
fprintf('Elektroda %d i elektroda %d imaju najlosiju povezanost (MSC :%f)\r', MinMSC(1),MinMSC(2),MinMSC(3)); 

%% Plotanje
%-----Histogram-----
figure();
im = imagesc(MSC_10Hz, [0 1]);
colorbar;
grid on;


%-----Plotanje pozicija elektroda (polarni)-----
figure()
for i = 1:16
    polarplot(deg2rad(elektroda{i}{2}),elektroda{i}{3}, 'b-o')
    text(deg2rad(elektroda{i}{2}),elektroda{i}{3}, [' ', char(elektroda{i}{4})]);
    hold on;
end

%-----Plotanje vr. razlike izmedju 2 signala sa najvecim MSC-om-----
figure('Name', 'MaxMSC')
plot(time, data_filt(MaxMSC(1),:))
hold on
plot(time, data_filt(MaxMSC(2),:))
grid on
title('Par elektroda s najvecim MSC-om')
xlim([200 210])

%-----Plotanje vr. razlike izmedju 2 signala sa najmanjim MSC-om-----
figure('Name', 'MinMSC')
plot(time, data_filt(MinMSC(1),:))
hold on
plot(time, data_filt(MinMSC(2),:))
grid on
title('Par elektroda s najmanjim MSC-om' )
xlim([200 210])
