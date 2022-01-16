clc, clear all, close all

%% import data: ICU Physiobank (8-bit)
Fs=125; 
PPG_3000060m=load("3000060m.mat"); PPG_3000060m=PPG_3000060m.val; 
PPG_3000063m=load("3000063m.mat"); PPG_3000063m=PPG_3000063m.val;
PPG_3000125m=load("3000125m.mat"); PPG_3000125m=PPG_3000125m.val;
PPG_3000142m=load("3000142m.mat"); PPG_3000142m=PPG_3000142m.val;
PPG_3000154m=load("3000154m.mat"); PPG_3000154m=PPG_3000154m.val;
PPG_3000397m=load("3000397m.mat"); PPG_3000397m=PPG_3000397m.val;
PPG_3000435m=load("3000435m.mat"); PPG_3000435m=PPG_3000435m.val;
PPG_3000701m=load("3000701m.mat"); PPG_3000701m=PPG_3000701m.val;
PPG_3000716m=load("3000716m.mat"); PPG_3000716m=PPG_3000716m.val;
PPG_3000912m=load("3000912m.mat"); PPG_3000912m=PPG_3000912m.val;
PPG_3001055m=load("3001055m.mat"); PPG_3001055m=PPG_3001055m.val;
PPG_3001133m=load("3001133m.mat"); PPG_3001133m=PPG_3001133m.val;

%% 
%create table
Patience=['PPG_3000060m';'PPG_3000063m'; 'PPG_3000125m'; 'PPG_3000142m';...
    'PPG_3000154m'; 'PPG_3000397m'; 'PPG_3000435m'; 'PPG_3000701m'; ...
    'PPG_3000716m'; 'PPG_3000912m'; 'PPG_3001055m'; 'PPG_3001133m'];
ICU_PhysioBank=[PPG_3000060m; PPG_3000063m; PPG_3000125m; PPG_3000142m; ...
    PPG_3000154m; PPG_3000397m; PPG_3000435m; PPG_3000701m; ...
    PPG_3000716m; PPG_3000912m; PPG_3001055m; PPG_3001133m];
PPG_table=table(Patience, ICU_PhysioBank); %create table
PPG=table2struct(PPG_table); %create struct
%%
N_Patiences=size(PPG,1);
N_Samples=size(PPG(1).ICU_PhysioBank,2);
time=(0:1/Fs:N_Samples/Fs-1/Fs)*10^3; %time in ms

% clear the single PPGs
clear PPG_3000060m PPG_3000063m PPG_3000125m PPG_3000142m...
    PPG_3000154m PPG_3000397m PPG_3000435m PPG_3000701m ...
    PPG_3000716m PPG_3000912m PPG_3001055m PPG_3001133m Patience...
    ICU_PhysioBank PPG_table

%visualization
figure()
for nPatience=1:N_Patiences
    subplot(3,4,nPatience), plot(time, PPG(nPatience).ICU_PhysioBank);
    xlabel('Time [ms]'), ylabel('PPG value'), title(nPatience);
end

%% Stage 1 - Clipping Bottom & Top
for nPatience=1:N_Patiences
    for nSample=1:N_Samples
        if PPG(nPatience).ICU_PhysioBank(nSample)<=0 || PPG(nPatience).ICU_PhysioBank(nSample)>=255
            PPG(nPatience).disturbed_indx(nSample)=1;
        end
    end
    if size(PPG(nPatience).ICU_PhysioBank,2)>=size(PPG(nPatience).disturbed_indx,2)
            PPG(nPatience).disturbed_indx=cat(2, PPG(nPatience).disturbed_indx(1:end), zeros(1, 1e6-size(PPG(nPatience).disturbed_indx,2)));
    end
end

%% Stage 2 & 3 - Filtering
%low pass
Fc_lp = 15;
Fc_hp=0.01;
[b_lp,a_lp] = butter(4,Fc_lp/(Fs/2),'low'); % butterworth filter, LP
[b_hp,a_hp] = butter(4,Fc_hp/(Fs/2),'high'); % butterworth filter, HP

%all signals
figure(),
for nPatience=1:N_Patiences
    PPG_lp=filtfilt(b_lp,a_lp, PPG(nPatience).ICU_PhysioBank); %LP
    PPG(nPatience).Filtered=filtfilt(b_hp,a_hp,PPG_lp); %LP+HP
    subplot(3,4,nPatience), plot(time, PPG(nPatience).Filtered);
    xlabel('Time [ms]'), ylabel('PPG value'), title(nPatience);
end

%comparison pre/post filtering
p=1;
figure(), 
plot(time, PPG(p).ICU_PhysioBank); hold on, 
plot(time, PPG(p).Filtered); hold off,
xlabel('Time [ms]'), ylabel('PPG value'), title('BP Filtering Effect');
legend('Original', 'Filtered')

%% cut the beginning of the signal
% p=12;
% figure(), plot(PPG(p).Filtered);
% xlabel('Sample'), ylabel('PPG value'), title(p)
figure(),
for nPatience=1:N_Patiences
    PPG(nPatience).Filtered_cut=PPG(nPatience).Filtered(3.8*10^4:end);
    if nPatience==1
        N_Samples_cut=size(PPG(nPatience).Filtered_cut,2);
        time_cut=(0:1/Fs:N_Samples_cut/Fs-1/Fs)*10^3; %time in ms
    end
     subplot(3,4,nPatience), plot(time_cut, PPG(nPatience).Filtered_cut);
     xlabel('Time [ms]'), ylabel('PPG value'), title(nPatience);
end

%% Stage 4: Peak-Valleys detection
%one signal only
y=PPG(6).Filtered_cut;
RingBuffer_signal=zeros(1, 4.8*Fs);
RingBuffer_annotations=zeros(1, 4.8*Fs);
RingBuffer_size=size(RingBuffer_signal,2);

figure(), plot(y(1:RingBuffer_size));
[X, Y]=ginput(3); %%selezionare minimi (PWB, PWE) e massimo (PWSP)
X=round(X); Y=round(Y);

PWD_pre=X(2)-X(1); %in samples
base=Y(1); % a mio parere si calcola sempre apex-valley(1) 
apex=Y(3);
%% Adaptive Threshold: MA filter on the last valid PWD (75%)
%initialization
B = (1/round(0.75*PWD_pre))*ones(round(0.75*PWD_pre),1);
AdThr = filter(B,1,y(X(1):X(2)+4.8*Fs)); %questo è un primo prototipo di AdThr, vorrei in realtà fosse lungo solo 4.8s !!!
delay=floor((round(0.75*PWD_pre)-1)/2); %delay?
first_pulse_length=X(2)-X(1)+1; %+1 little correction

%%
%search for absolute max and min in the next pulse
%peaks
k=1;
pos_start(k) = X(2);
%pos_start dovrà essere aggiornato ciclicamente, per ogni finestra
%e' la valle successiva al pulse appena osservato
%POS_START è SEMPRE uguale all'ultima valle che considero, ossia l'inizio
%del ring buffer
[peak_candidate, loc_peak]=findpeaks(y(pos_start(k):pos_start(k)+4.8*Fs)); %locs mi da il punto dove c'è il peak, a partire da dove conto
%loc_peak_candidate = pos_start(k)+loc_peak; %findpeaks legge solo a partire dal valore pos_start, in questo modo risolvo la traslazione

loc_peak_AdThr=zeros(1,length(loc_peak)); %preallocating space for speed
j=1;
for i=1:length(loc_peak)
    if peak_candidate(i)>AdThr(first_pulse_length+loc_peak(i)) %se valore(primo peak) > valore(Adthr) : peak(1)>Adthr(posizione_peak(1))
        loc_peak_AdThr(j)=loc_peak(i); %attenzione! usano due indici diversi, perchè il numero di peak_validi sarà <= peak_trovati
        j=j+1;
    end
end
loc_peak_AdThr = loc_peak_AdThr(1:j-1)
%%
[valley_candidate, loc_valley]=findpeaks(-y(pos_start(k):pos_start(k)+4.8*Fs));
%loc_valley_candidate= pos_start(k)+loc_valley;

%This could be needed just in the FIRST ITERATION->CHECK
loc_valley=[0 loc_valley];
valley_candidate=[y(pos_start(k)) valley_candidate];
% 
loc_valley_AdThr=zeros(1,length(loc_valley));

j=1;
for i=1:length(loc_valley)
    if valley_candidate(i)>-AdThr(first_pulse_length+loc_valley(i)) %se valore(primo peak) > valore(Adthr) : peak(1)>Adthr(posizione_peak(1))
        loc_valley_AdThr(j)=loc_valley(i);
        j = j+1;
    end
end

loc_valley_AdThr = loc_valley_AdThr(1:j-1);
loc_valley_AdThr = [0 loc_valley_AdThr]; %ZERO?
%%
figure(),
plot(y(pos_start(k):pos_start(k)+4.8*Fs)), hold on,
plot(AdThr(first_pulse_length:end), '--r'),
%plot(first_pulse_length+loc_peak_AdThr-1, y(pos_start(k)+loc_peak_AdThr-1), '*m'), 
plot(loc_valley_AdThr, y(pos_start(k)+loc_valley_AdThr), '*g'),
%
plot(loc_peak, y(pos_start(k)+loc_peak), '*r')

%////////////////////////////////////////////////////////////%////////////////////////////////////////////////////////////%////////////////////////////////////////////////////////////
%%
%TEMPORARY SOLUTION: we make this fix just for dicrotic notches that appear
%under AdThreshold curve. We should also do this for the ones appearing
%above.

%Afterwards, we merge these outcomes with what we found out with AdThr
%comparison
proto_length=min(length(loc_valley),length(loc_peak));
PWA_proto = zeros(1,proto_length);
for i=1:proto_length % I could just grab less peak-valley couplesssss, e.g. 5
    PWA_proto(i)=y(pos_start(k)+loc_peak(i))-y(pos_start(k)+loc_valley(i));
end 

%so I start from the second element
j=2;
m=1;

PWA=zeros(1,length(PWA_proto)-1);
loc_valley_PWA = zeros(1,length(PWA_proto));
loc_peak_PWA = zeros(1,length(PWA_proto));
loc_peak_dicro = zeros(1,length(PWA_proto));
loc_valley_dicro = zeros(1,length(PWA_proto));

%first one is always a pulse,as it comes from previous ringbuffer's
%considerations

%/////////////////////////
PWA(1)=PWA_proto(1);
loc_valley_PWA(1)=loc_valley(1); %first valley is always starting point of ringbuffer
loc_peak_PWA(1)=loc_peak(1); % ATTENZIONE QUESTO TEMPORARY FIX, WILL IT WORK LATER?

for i=2:length(PWA_proto)
    if (PWA_proto(i)*3>PWA_proto(i-1)) %this discriminates btw real pulses and dicrotic peaks
        PWA(j)=PWA_proto(i);
        loc_valley_PWA(j)=loc_valley(i);
        loc_peak_PWA(j)=loc_peak(i);
        j=j+1;
    else
        loc_peak_dicro(m)=loc_peak(i);
        loc_valley_dicro(m)=loc_valley(i);
        m=m+1;
    end
end
%/////////////////////////
%I have initialized these being generous with size, now I cut them short
PWA=PWA(1:j-1);
loc_valley_PWA = loc_valley_PWA(1:j-1);
loc_peak_PWA = loc_peak_PWA(1:j-1);
loc_peak_dicro = loc_peak_dicro(1:m-1);
loc_valley_dicro = loc_valley_dicro(1:m-1);
%%

%What I want to have is a CLEAR discrimination btw PULSE peak/valley elements and
%segments, and DICROTIC peak/valley.
%What I'm interested into is the X local (1-600) location of  both peaks (and
%later I will look into valleys)

%So PEAKS have to double check AdThr and PWA3
%DICROs have already been looked into

loc_peak_merge = intersect(loc_peak_AdThr, loc_peak_PWA);
loc_valley_merge = intersect(loc_valley_AdThr, loc_valley_PWA);
%TEMPORARY
loc_valley_merge = [loc_valley_merge]; %ZERO?

figure(),
title('POST_MERGE');
plot(y(pos_start(k):pos_start(k)+4.8*Fs)), hold on,
plot(AdThr(first_pulse_length:end), '--r'), 
plot(loc_peak_merge, y(pos_start(k)+loc_peak_merge), '*m'), 
plot(loc_valley_merge, y(pos_start(k)+loc_valley_merge), '*g'),
%
plot(loc_peak_dicro, y(pos_start(k)+loc_peak_dicro), '*r')


%% Checks
artifact=[]
while(1)
    %Check3
    PWA_left=(y(pos_start(k)+loc_peak_merge(1))-y(pos_start(k)+loc_valley_merge(1)));
    PWA_right=(y(pos_start(k)+loc_peak_merge(1))-y(pos_start(k)+loc_valley_merge(2)));
    if(PWA_left<=2*mean(abs(diff(y(pos_start(k):pos_start(k)+4.8*Fs))))) %Check3 TRUE
        z=3;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check3 FALSE
    
    %Check4
    min_time = 0.08*Fs; %time in samples
    max_time = 0.49*Fs;
    PWRT=loc_peak_merge(1)-loc_valley_merge(1);
    if(PWRT<min_time || PWRT>max_time) %Check4 TRUE
        z=4;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check4 FALSE
    
    %Check5
    SisTime=(loc_peak_merge(1)-loc_valley_merge(1));
    DiasTime=(loc_valley_merge(2)-loc_peak_merge(1));
    PWSDratio=SisTime/DiasTime;
    if(PWSDratio>1.1) %Check5 TRUE
        z=5;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check5 FALSE
    
    %Check6
    min_time = 0.27*Fs; %time in samples
    max_time = 2.4*Fs;
    PWD=loc_valley_merge(2)-loc_valley_merge(1);
    
    if(PWD<min_time || PWD>max_time) %Check6 TRUE
        z=6;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check6 FALSE
    
    %Check7
    dc=0;
    for i=1:length(loc_peak_dicro)
        if(loc_peak_dicro(i)>loc_peak_merge(1) && loc_peak_dicro(i) <loc_valley_merge(2))
            dc=dc+1;
        end
    end
    
    if(dc>2) %Check6 TRUE % Shall this be >2 or >=2???
        z=7;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check7 FALSE
    
    %Check8
    for i=1:PWRT
        if(y(pos_star(k)+loc_valley_merge(1)+1)<=y(pos_star(k)+loc_valley_merge(1)))
            z=8;
            artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
            pos_start(k)=pos_start(k)+loc_valley_merge(2);
            break
        end
    end
        
    %Check8 FALSE  
    
    %Check9
    f=1;
    for i=1:length(loc_valley_dicro)
        if(loc_valley_dicro(i)>loc_peak_merge(1) && loc_valley_dicro(i) <loc_valley_merge(2)) %I save the valleys (aka lower points) btw peak1 and valley2
            loc_valley_dicro_fp(f)=loc_valley_dicro(i);
            f=f+1;
        end
    end
    
    for i=1:f-1
        if(y(loc_valley_dicro_fp)<y(loc_valley_merge(2)
            z=9;
            artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
            pos_start(k)=pos_start(k)+loc_valley_merge(2);
            break
        end
    end
    
    %Check9 FALSE
    
    %Check10
    PWA_LR=PWA_left/PWA_right;
    PWA_RL=1/PWA_LR;
    if(PWA_LR <0.4 || PWA_RL <0.4)
        z=10;
        artifact=[artifact pos_start(k)+loc_valley_merge(1) pos_start(k)+loc_valley_merge(2)-1];
        pos_start(k)=pos_start(k)+loc_valley_merge(2);
        break
    end
    
    %Check10 FALSE
    
    %Check11
    
    
    
    
    
    
    
    %
    pos_start(k)=pos_start(k)+loc_valley_merge(2);
    break
end

