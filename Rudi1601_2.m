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
y=PPG(1).Filtered_cut;
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
%
% figure(),
% plot(y(1:PWD_pre-1)), hold on, plot(AdThr(1:PWD_pre-1));
% 
% figure(), plot(y(X(1):X(2)-1)), hold on, 
% plot(0:length(AdThr(X(1):X(2)))-delay-1, AdThr(X(1):X(2)-delay), '--r'), hold off,
% xlabel('Sample'), title('MA of signal'), legend('Signal', 'MA-delay');

% y_filtered=filter(B,1,y);
% figure(), plot (y), hold on, 
% plot(0:length(y_filtered)-delay-1, y_filtered(1:end-delay), '--r'), hold off;
% xlabel('Sample'), title('MA of signal'), legend('Signal', 'MA-delay');
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

loc_valley_AdThr=zeros(1,length(loc_valley));

j=1;
for i=1:length(loc_valley)
    if valley_candidate(i)>-AdThr(first_pulse_length+loc_valley(i)) %se valore(primo peak) > valore(Adthr) : peak(1)>Adthr(posizione_peak(1))
        loc_valley_AdThr(j)=loc_valley(i);
        j = j+1;
    end
end

loc_valley_AdThr = loc_valley_AdThr(1:j-1);
loc_valley_AdThr = [0 loc_valley_AdThr];
%%
figure(),
plot(y(X(1): pos_start(k)+4.8*Fs)), hold on,
plot(AdThr, '--r'), 
%plot(first_pulse_length+loc_peak_AdThr-1, y(pos_start(k)+loc_peak_AdThr-1), '*m'), 
plot(first_pulse_length+loc_valley_AdThr-1, y(pos_start(k)+loc_valley_AdThr-1), '*g'),
%
plot(first_pulse_length+loc_peak-1, y(pos_start(k)+loc_peak-1), '*r')

%////////////////////////////////////////////////////////////%////////////////////////////////////////////////////////////%////////////////////////////////////////////////////////////
%%
%TEMPORARY SOLUTION: we make this fix just for dicrotic notches that appear
%under AdThreshold curve. We should also do this for the ones appearing
%above.

%Afterwards, we merge these outcomes with what we found out with AdThr
%comparison

PWA_proto = zeros(1,length(loc_peak));
for i=1:length(loc_peak) % I could just grab less peak-valley couplesssss, e.g. 5
    PWA_proto(i)=y(pos_start(k)+loc_peak(i))-y(pos_start(k)+loc_valley_AdThr(i));
end 

%so I start from the second element
j=2;
m=1;

PWA=zeros(1,length(PWA_proto)-1);
loc_valley_PWA = zeros(1,length(PWA_proto)-1);
loc_peak_PWA = zeros(1,length(PWA_proto)-1);
loc_peak_dicro = zeros(1,length(PWA_proto)-1);


%first one is always a pulse,as it comes from previous ringbuffer's
%considerations

PWA(1)=PWA_proto(1);
loc_valley_PWA(1)=loc_valley_AdThr(1);
loc_peak_PWA(1)=loc_peak(1); % ATTENZIONE QUESTO TEMPORARY FIX, WILL IT WORK LATER?

for i=2:length(PWA_proto)
    if (PWA_proto(i)*3>PWA_proto(i-1)) %this discriminates btw real pulses and dictotic peaks
        PWA(j)=PWA_proto(i);
        loc_valley_PWA(j)=loc_valley_AdThr(i);
        loc_peak_PWA(j)=loc_peak(i);
        j=j+1;
    else
        loc_peak_dicro(m)=loc_peak(i);
        m=m+1;
    end
end

%I have initialized these being generous with size, now I cut them short
PWA=PWA(1:j-1);
loc_valley_PWA = loc_valley_PWA(1:j-1);
loc_peak_PWA = loc_peak_PWA(1:j-1);
loc_peak_dicro = loc_peak_dicro(1:m-1);

%%

%What I want to have is a CLEAR discrimination btw PULSE peak/valley elements and
%segments, and DICROTIC peak/valley.
%What I'm interested into is the X local (1-600) location of  both peaks (and
%later I will look into valleys)

%So PEAKS have to double check AdThr and PWA3
%DICROs have already been looked into

loc_peak_merge = intersect(loc_peak_AdThr, loc_peak_PWA);

figure(),
plot(y(X(1): pos_start(k)+4.8*Fs)), hold on,
plot(AdThr, '--r'), 
plot(first_pulse_length+loc_peak_merge-1, y(pos_start(k)+loc_peak_merge-1), '*m'), 
plot(first_pulse_length+loc_valley_AdThr-1, y(pos_start(k)+loc_valley_AdThr-1), '*g'),
%
plot(first_pulse_length+loc_peak_dicro-1, y(pos_start(k)+loc_peak_dicro-1), '*r')


%% Check 3
% There is an incomprension: PWA_pre is considered a reference in our
% calculation, but in this case we're doing condition on it, not on the
% actual signal...

%That said, I'm not sure whether I'm considering N-1 as a group of
%parameters of the previous pulse wave I'm looking into, or as the first
%pulse wave of the actual Ring Buffer.

%Done these considerations, I think I'll consider it as the first pulse
%wave of the buffer, while N-0 is a pulse wave being recorded in THIS
%moment, so not yet subject on the filtering.
%All the conditions will be applied on N-1, so I think that is coherent.
%mean(abs(diff(RingBufferSignal)

Check_3_value=mean(abs(y(pos_start:pos_start+4.8*Fs))); % i have to add diff

if (y(loc_peak_AdThr(1))-y(loc_valley_AdThr(1)))< (2*Check_3_value)
    %se minore:artifact: sets, 2 at a time, the start & the end of an
    %artifact area; there surely is a more clever way to mark it, will make
    %it coherent later
    artifact=[artifact loc_valley_AdThr(1) loc_valley_AdThr(2)];
end
%break?
%% Check 4
%Time between valley(1) and peak(1) shall be more than 0.08s and less 0.49s

min_time = 0.08*Fs; %time in samples
max_time = 0.49*Fs;
if(loc_peak_AdThr(1)-loc_valley_AdThr(1) < min_time || loc_peak_AdThr(1)-loc_valley_AdThr(1) > max_time)
    artifact=[artifact loc_valley_AdThr(1) loc_valley_AdThr(2)];
end

%% Check 5
%SistolicTime/DiastolicTime not bigger than 1.1

SisTime=loc_peak_AdThr(1)-loc_valley_AdThr(1);
DiasTime=loc_valley_AdThr(2)-loc_peak_AdThr(1);

if(SisTime/DiasTime >1.1)
    artifact=[artifact loc_valley_AdThr(1) loc_valley_AdThr(2)];
end
