%%%%%%%%%%%to run it
% clear all
% close all
% directory1 = ['D:\TDT_mmn paradigms\Intan\RHD_MATLAB_functions\'];%folder to save the files to
% directory2 = ['D:\TDT_mmn paradigms\Intan\RHD_MATLAB_functions\']; % this is where the TDT event file is 
% Chs1 = [2:17];%Amp A chs:2-32
% Chs2 = [34:49];%AmpB chs: 34-64
% epoch_tframe = [-25 200];
% read_Intan_RHD2000_file
% 
% 
% filenamesout=[filename(1:end-10) '@ep' '.mat'];
%  te1=board_dig_in_data(1,:);
%  AA=amplifier_data;
% adrate = frequency_parameters.amplifier_sample_rate;
% clearvars -except amplifier_channels amplifier_data board_dig_in_data frequency_parameters filename directory1 directory2 Chs1 Chs2 epoch_tframe te1 AA adrate filenamesout
% tic
% [eegm1,eege1,eegm2,eege2,time] =ImportIntanMMN2electrodes(directory1,directory2,Chs1,Chs2,epoch_tframe,te1,AA,adrate,filenamesout);
% toc

function  [eegm1,eege1,eegm2,eege2,time] =ImportIntanMMN2electrodes(directory1,directory2,Chs1,Chs2,epoch_tframe,te1,AA,adrate,filenamesout);
     newadrate = 1000; 
 filenames0   = {};
 if isempty(filenames0)
    filenames=[];
    f=dir( [ directory2 '*.csv']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
 end
  %to import the data in the .csv (event files)
 x=importdata([ directory1 filenames{1} ]); %data re standards and Devs

 %get trigger info from Digital input
rising_edges = find(diff(te1 > .1) == 1); %this is the times of the trigger onsets
%%get out the triggers for MMN
D=find(x.data(:,13) ==2);%these are the locations of the deviant triggers
S=find(x.data(:,13) ==0);%these are the locations of the Standard triggers
ss=D-1;%these are the locations of the Standards right before the devinats
if x.data(:,10)==0
    ST="tone";
elseif x.data(:,10)==1
    ST="AM";
else x.data(:,10)==2
    ST="FM";
end

%%% Resample trigger data
for i1=1:length(rising_edges)
    trig01(i1)    = round(rising_edges(i1)./(adrate/newadrate));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR Electrode on AMP A
%get out raw data
% AA=amplifier_data;
X=AA(Chs1(1):Chs1(end),:); 
 %%%%%to filter the raw data into LFP and MUA
             trigtypesel = [];
             newadrate = 1000; 
             
             filtere = [0.5 300];
             filteru = [300 5000];
             filtertype=1;
             fsize = 6;
             xlabelres=5;
             time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);
%               adrate = frequency_parameters.amplifier_sample_rate;           
        cnte = zeros(size(X));
        cntm = zeros(size(X));
        cntu = zeros(size(X));

        % filtering to get field
        if filtere(1)==0
            n = 2;
            Wn = filtere(2)/(adrate/2);
            [b,a] = butter(n,Wn,'low');
            for i=1:size(X,1)
                cnte(i,:)=filtfilt(b,a,X(i,:));
            end
        else
            n = 2;
            Wn = filtere/(adrate/2);
            [b,a] = butter(n,Wn);
            for i=1:size(X,1)
                cnte(i,:)=filtfilt(b,a,X(i,:));
            end
        end
        %%%% filtering to get the unit/MUA
        if filteru(2)==0
            n = 2;
            Wn = filteru(1)/(adrate/2);
            [b,a] = butter(n,Wn,'high');
            for i=1:size(X,1)
                cntu(i,:)=filtfilt(b,a,X(i,:));
            end
        else
            n = 2;
            Wn = filteru/(adrate/2);
            [b,a] = butter(n,Wn);
            for i=1:size(X,1)
                cntu(i,:)=filtfilt(b,a,X(i,:));
            end
        end
          % rectifying unit
        for i=1:size(cntu,1)
            cntm(i,:)=abs(hilbert(cntu(i,:)));
        end


            if adrate ~= newadrate
    
                % filtering before downsample
                nyq = newadrate/2.01;
    
                if filtertype == 1
    
    
                    n = 2;
                    Wn = nyq/(adrate/2);
                    [b,a] = butter(n,Wn,'low');
                    for i=1:size(cnte,1)
                        cnte(i,:)=filtfilt(b,a,cnte(i,:));
    
    
                    end
                    for i=1:size(cntu,1)
    
                        cntm(i,:)=filtfilt(b,a,cntm(i,:));
    
                        cntu(i,:)=filtfilt(b,a,cntu(i,:));
    
                    end
    
                else
    
                    bpFilt = designfilt('lowpassfir','FilterOrder',200, ...
                        'CutoffFrequency',nyq, ...
                        'SampleRate',adrate);
                    for i=1:size(cnte,1)
                        cnte(i,:)=filtfilt(bpFilt,cnte(i,:));
    
    
                    end
                    for i=1:size(cntu,1)
    
                        cntm(i,:)=filtfilt(bpFilt,cntm(i,:));
    
                        cntu(i,:)=filtfilt(bpFilt,cntu(i,:));
    
                    end
                end
    
    
                % downsample or resmaple
                downsampleby    = adrate/newadrate;
                if downsampleby-round(downsampleby)==0
                    cnte            = downsample(cnte',downsampleby)';
    
                    cntm            = downsample(cntm',downsampleby)';
                    cntu            = downsample(cntu',downsampleby)';
    
                    try
                        cnt_arej        = downsample(craw.arej',downsampleby)';
                    catch
                        cnt_arej        = [];
                    end
                else
                    cnte            = resample(cnte',newadrate,adrate)';
    
                    cntm            = resample(cntm',newadrate,adrate)';
                    cntu            = resample(cntu',newadrate,adrate)';
    
                    try
                        cnt_arej        = resample(craw.arej',newadrate,adrate)';
                    catch
                        cnt_arej        = [];
                    end
                end
    
            else
    
                try
                    size(cnt_arej);
                catch
                    cnt_arej        = [];
                end
            end
             %%%Epoching the data 
      x1 = round(epoch_tframe(1)*(newadrate/1000));
      x2 = round(epoch_tframe(2)*(newadrate/1000));
      for i=1:length(D) %this is for the Dev stimuli
          eegDm(:,i,:)=cntm(:,trig01(D(i))+x1:trig01(D(i))+x2);%epoching the mua
          eegDe(:,i,:)=cnte(:,trig01(D(i))+x1:trig01(D(i))+x2);%epoching the LFP
      end
        for i=1:length(ss)
          eegSm(:,i,:)=cntm(:,trig01(ss(i))+x1:trig01(ss(i))+x2);%epoching the mua
          eegSe(:,i,:)=cnte(:,trig01(ss(i))+x1:trig01(ss(i))+x2);%epoching the LFP
      end
          % baseline before image
          clear x1 x2
          x1 = find(time<=epoch_tframe(1), 1, 'last' );
          x2 = find(time<=-5, 1, 'last' );
    
          for i=1:size(eegSm,1)
              for iii=1:size(eegSm,2)
                  eegSm(i,iii,:)=squeeze(eegSm(i,iii,:))-squeeze(mean(eegSm(i,iii,x1:x2),3));
                  eegDm(i,iii,:)=squeeze(eegDm(i,iii,:))-squeeze(mean(eegDm(i,iii,x1:x2),3));
                  eegDe(i,iii,:)=squeeze(eegDe(i,iii,:))-squeeze(mean(eegDe(i,iii,x1:x2),3));
                  eegSe(i,iii,:)=squeeze(eegSe(i,iii,:))-squeeze(mean(eegSe(i,iii,x1:x2),3));
              end
          end
    
          avgDe=squeeze(mean(eegDe(:,:,:),2));
          avgDm=squeeze(mean(eegDm(:,:,:),2));
          avgSe=squeeze(mean(eegSe(:,:,:),2));
          avgSm=squeeze(mean(eegSm(:,:,:),2));


          %%To make overlay responses for stnad & dev for each ch
           curfig = figure;
                    set(curfig,'position',[100   100   1500  1500],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                    subplot(16,2,1)
                    plot(time,avgSe(1,:))
                    hold on
                    plot(time,avgDe(1,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch1, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,2)
                    plot(time,avgSm(1,:))
                    hold on
                    plot(time,avgDm(1,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch1, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,3)
                    plot(time,avgSe(2,:))
                    hold on
                    plot(time,avgDe(2,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch2, B=stand, R=Dev' ] )
                    subplot(16,2,4)
                    plot(time,avgSm(2,:))
                    hold on
                    plot(time,avgDm(2,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch2, B=stand, R=Dev' ] )
                      subplot(16,2,5)
                    plot(time,avgSe(3,:))
                    hold on
                    plot(time,avgDe(3,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch3, B=stand, R=Dev' ] )
                    subplot(16,2,6)
                    plot(time,avgSm(3,:))
                    hold on
                    plot(time,avgDm(3,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch3, B=stand, R=Dev' ] )
                      subplot(16,2,7)
                    plot(time,avgSe(4,:))
                    hold on
                    plot(time,avgDe(4,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch4, B=stand, R=Dev' ] )
                    subplot(16,2,8)
                    plot(time,avgSm(4,:))
                    hold on
                    plot(time,avgDm(4,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch4, B=stand, R=Dev' ] )
                      subplot(16,2,9)
                    plot(time,avgSe(5,:))
                    hold on
                    plot(time,avgDe(5,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch5, B=stand, R=Dev' ] )
                    subplot(16,2,10)
                    plot(time,avgSm(5,:))
                    hold on
                    plot(time,avgDm(5,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch5, B=stand, R=Dev' ] )
                    subplot(16,2,11)
                    plot(time,avgSe(6,:))
                    hold on
                    plot(time,avgDe(6,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch6, B=stand, R=Dev' ] )
                    subplot(16,2,12)
                    plot(time,avgSm(6,:))
                    hold on
                    plot(time,avgDm(6,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch6, B=stand, R=Dev' ] )
                    subplot(16,2,13)
                    plot(time,avgSe(7,:))
                    hold on
                    plot(time,avgDe(7,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch7, B=stand, R=Dev' ] )
                    subplot(16,2,14)
                    plot(time,avgSm(7,:))
                    hold on
                    plot(time,avgDm(7,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch7, B=stand, R=Dev' ] )
                    subplot(16,2,15)
                    plot(time,avgSe(8,:))
                    hold on
                    plot(time,avgDe(8,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch8. B=stand, R=Dev' ] )
                    subplot(16,2,16)
                    plot(time,avgSm(8,:))
                    hold on
                    plot(time,avgDm(8,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch8, B=stand, R=Dev' ] )

                    subplot(16,2,17)
                    plot(time,avgSe(9,:))
                    hold on
                    plot(time,avgDe(9,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch9, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,18)
                    plot(time,avgSm(9,:))
                    hold on
                    plot(time,avgDm(9,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch9, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,19)
                    plot(time,avgSe(10,:))
                    hold on
                    plot(time,avgDe(10,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch10, B=stand, R=Dev' ] )
                    subplot(16,2,20)
                    plot(time,avgSm(10,:))
                    hold on
                    plot(time,avgDm(10,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch10, B=stand, R=Dev' ] )
                      subplot(16,2,21)
                    plot(time,avgSe(11,:))
                    hold on
                    plot(time,avgDe(11,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch11, B=stand, R=Dev' ] )
                    subplot(16,2,22)
                    plot(time,avgSm(11,:))
                    hold on
                    plot(time,avgDm(11,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch11, B=stand, R=Dev' ] )
                      subplot(16,2,23)
                    plot(time,avgSe(12,:))
                    hold on
                    plot(time,avgDe(12,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch12, B=stand, R=Dev' ] )
                    subplot(16,2,24)
                    plot(time,avgSm(12,:))
                    hold on
                    plot(time,avgDm(12,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch12, B=stand, R=Dev' ] )
                      subplot(16,2,25)
                    plot(time,avgSe(13,:))
                    hold on
                    plot(time,avgDe(13,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch13, B=stand, R=Dev' ] )
                    subplot(16,2,26)
                    plot(time,avgSm(13,:))
                    hold on
                    plot(time,avgDm(13,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch13, B=stand, R=Dev' ] )
                    subplot(16,2,27)
                    plot(time,avgSe(14,:))
                    hold on
                    plot(time,avgDe(14,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch14, B=stand, R=Dev' ] )
                    subplot(16,2,28)
                    plot(time,avgSm(14,:))
                    hold on
                    plot(time,avgDm(14,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch14, B=stand, R=Dev' ] )
                    subplot(16,2,29)
                    plot(time,avgSe(15,:))
                    hold on
                    plot(time,avgDe(15,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch15, B=stand, R=Dev' ] )
                    subplot(16,2,30)
                    plot(time,avgSm(15,:))
                    hold on
                    plot(time,avgDm(15,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch15, B=stand, R=Dev' ] )
                    subplot(16,2,31)
                    plot(time,avgSe(16,:))
                    hold on
                    plot(time,avgDe(16,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch16. B=stand, R=Dev' ] )
                    subplot(16,2,32)
                    plot(time,avgSm(16,:))
                    hold on
                    plot(time,avgDm(16,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch16, B=stand, R=Dev' ] )

                        axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,['1-' filenamesout  '-'  ', ' num2str(length(trig01)) ' sweeps'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                     print( '-djpeg','-r300',[directory1 '1-LFP__MUA_Stnd_Dev Stim Type=' num2str(ST) '_' filenamesout  '.jpg'])

avgDe1=avgDe;
avgDm1=avgDm;
avgSe1=avgSe;
avgSm1=avgSm;
XA=X;
clear avgSm avgSe X avgDm avgDm eegSm eegDm eegDe eegSe cntm cnte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FOR Electrode on AMP B
%get out raw data
X=AA(Chs2(1):Chs2(end),:); 
%%%%%to filter the raw data into LFP and MUA
             trigtypesel = [];
             newadrate = 1000; 
             
             filtere = [0.5 300];
             filteru = [300 5000];
             filtertype=1;
             fsize = 6;
             xlabelres=5;
             time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);
%               adrate = frequency_parameters.amplifier_sample_rate;           
        cnte = zeros(size(X));
        cntm = zeros(size(X));
        cntu = zeros(size(X));

        % filtering to get field
        if filtere(1)==0
            n = 2;
            Wn = filtere(2)/(adrate/2);
            [b,a] = butter(n,Wn,'low');
            for i=1:size(X,1)
                cnte(i,:)=filtfilt(b,a,X(i,:));
            end
        else
            n = 2;
            Wn = filtere/(adrate/2);
            [b,a] = butter(n,Wn);
            for i=1:size(X,1)
                cnte(i,:)=filtfilt(b,a,X(i,:));
            end
        end
        %%%% filtering to get the unit/MUA
        if filteru(2)==0
            n = 2;
            Wn = filteru(1)/(adrate/2);
            [b,a] = butter(n,Wn,'high');
            for i=1:size(X,1)
                cntu(i,:)=filtfilt(b,a,X(i,:));
            end
        else
            n = 2;
            Wn = filteru/(adrate/2);
            [b,a] = butter(n,Wn);
            for i=1:size(X,1)
                cntu(i,:)=filtfilt(b,a,X(i,:));
            end
        end
          % rectifying unit
        for i=1:size(cntu,1)
            cntm(i,:)=abs(hilbert(cntu(i,:)));
        end


            if adrate ~= newadrate
    
                % filtering before downsample
                nyq = newadrate/2.01;
    
                if filtertype == 1
    
    
                    n = 2;
                    Wn = nyq/(adrate/2);
                    [b,a] = butter(n,Wn,'low');
                    for i=1:size(cnte,1)
                        cnte(i,:)=filtfilt(b,a,cnte(i,:));
    
    
                    end
                    for i=1:size(cntu,1)
    
                        cntm(i,:)=filtfilt(b,a,cntm(i,:));
    
                        cntu(i,:)=filtfilt(b,a,cntu(i,:));
    
                    end
    
                else
    
                    bpFilt = designfilt('lowpassfir','FilterOrder',200, ...
                        'CutoffFrequency',nyq, ...
                        'SampleRate',adrate);
                    for i=1:size(cnte,1)
                        cnte(i,:)=filtfilt(bpFilt,cnte(i,:));
    
    
                    end
                    for i=1:size(cntu,1)
    
                        cntm(i,:)=filtfilt(bpFilt,cntm(i,:));
    
                        cntu(i,:)=filtfilt(bpFilt,cntu(i,:));
    
                    end
                end
    
    
                % downsample or resmaple
                downsampleby    = adrate/newadrate;
                if downsampleby-round(downsampleby)==0
                    cnte            = downsample(cnte',downsampleby)';
    
                    cntm            = downsample(cntm',downsampleby)';
                    cntu            = downsample(cntu',downsampleby)';
    
                    try
                        cnt_arej        = downsample(craw.arej',downsampleby)';
                    catch
                        cnt_arej        = [];
                    end
                else
                    cnte            = resample(cnte',newadrate,adrate)';
    
                    cntm            = resample(cntm',newadrate,adrate)';
                    cntu            = resample(cntu',newadrate,adrate)';
    
                    try
                        cnt_arej        = resample(craw.arej',newadrate,adrate)';
                    catch
                        cnt_arej        = [];
                    end
                end
    
            else
    
                try
                    size(cnt_arej);
                catch
                    cnt_arej        = [];
                end
            end
             %%%Epoching the data 
      x1 = round(epoch_tframe(1)*(newadrate/1000));
      x2 = round(epoch_tframe(2)*(newadrate/1000));
      for i=1:length(D) %this is for the Dev stimuli
          eegDm(:,i,:)=cntm(:,trig01(D(i))+x1:trig01(D(i))+x2);%epoching the mua
          eegDe(:,i,:)=cnte(:,trig01(D(i))+x1:trig01(D(i))+x2);%epoching the LFP
      end
        for i=1:length(ss)
          eegSm(:,i,:)=cntm(:,trig01(ss(i))+x1:trig01(ss(i))+x2);%epoching the mua
          eegSe(:,i,:)=cnte(:,trig01(ss(i))+x1:trig01(ss(i))+x2);%epoching the LFP
      end
          % baseline before image
          clear x1 x2
          x1 = find(time<=epoch_tframe(1), 1, 'last' );
          x2 = find(time<=-5, 1, 'last' );
    
          for i=1:size(eegSm,1)
              for iii=1:size(eegSm,2)
                  eegSm(i,iii,:)=squeeze(eegSm(i,iii,:))-squeeze(mean(eegSm(i,iii,x1:x2),3));
                  eegDm(i,iii,:)=squeeze(eegDm(i,iii,:))-squeeze(mean(eegDm(i,iii,x1:x2),3));
                  eegDe(i,iii,:)=squeeze(eegDe(i,iii,:))-squeeze(mean(eegDe(i,iii,x1:x2),3));
                  eegSe(i,iii,:)=squeeze(eegSe(i,iii,:))-squeeze(mean(eegSe(i,iii,x1:x2),3));
              end
          end
    
          avgDe=squeeze(mean(eegDe(:,:,:),2));
          avgDm=squeeze(mean(eegDm(:,:,:),2));
          avgSe=squeeze(mean(eegSe(:,:,:),2));
          avgSm=squeeze(mean(eegSm(:,:,:),2));


          %%To make overlay responses for stnad & dev for each ch
           curfig = figure;
                    set(curfig,'position',[100   100   1500  1500],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                    subplot(16,2,1)
                    plot(time,avgSe(1,:))
                    hold on
                    plot(time,avgDe(1,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch1, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,2)
                    plot(time,avgSm(1,:))
                    hold on
                    plot(time,avgDm(1,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch1, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,3)
                    plot(time,avgSe(2,:))
                    hold on
                    plot(time,avgDe(2,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch2, B=stand, R=Dev' ] )
                    subplot(16,2,4)
                    plot(time,avgSm(2,:))
                    hold on
                    plot(time,avgDm(2,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch2, B=stand, R=Dev' ] )
                      subplot(16,2,5)
                    plot(time,avgSe(3,:))
                    hold on
                    plot(time,avgDe(3,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch3, B=stand, R=Dev' ] )
                    subplot(16,2,6)
                    plot(time,avgSm(3,:))
                    hold on
                    plot(time,avgDm(3,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch3, B=stand, R=Dev' ] )
                      subplot(16,2,7)
                    plot(time,avgSe(4,:))
                    hold on
                    plot(time,avgDe(4,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch4, B=stand, R=Dev' ] )
                    subplot(16,2,8)
                    plot(time,avgSm(4,:))
                    hold on
                    plot(time,avgDm(4,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch4, B=stand, R=Dev' ] )
                      subplot(16,2,9)
                    plot(time,avgSe(5,:))
                    hold on
                    plot(time,avgDe(5,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch5, B=stand, R=Dev' ] )
                    subplot(16,2,10)
                    plot(time,avgSm(5,:))
                    hold on
                    plot(time,avgDm(5,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch5, B=stand, R=Dev' ] )
                    subplot(16,2,11)
                    plot(time,avgSe(6,:))
                    hold on
                    plot(time,avgDe(6,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch6, B=stand, R=Dev' ] )
                    subplot(16,2,12)
                    plot(time,avgSm(6,:))
                    hold on
                    plot(time,avgDm(6,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch6, B=stand, R=Dev' ] )
                    subplot(16,2,13)
                    plot(time,avgSe(7,:))
                    hold on
                    plot(time,avgDe(7,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch7, B=stand, R=Dev' ] )
                    subplot(16,2,14)
                    plot(time,avgSm(7,:))
                    hold on
                    plot(time,avgDm(7,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch7, B=stand, R=Dev' ] )
                    subplot(16,2,15)
                    plot(time,avgSe(8,:))
                    hold on
                    plot(time,avgDe(8,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch8. B=stand, R=Dev' ] )
                    subplot(16,2,16)
                    plot(time,avgSm(8,:))
                    hold on
                    plot(time,avgDm(8,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch8, B=stand, R=Dev' ] )

                    subplot(16,2,17)
                    plot(time,avgSe(9,:))
                    hold on
                    plot(time,avgDe(9,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch9, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,18)
                    plot(time,avgSm(9,:))
                    hold on
                    plot(time,avgDm(9,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch9, B=stand, R=Dev, Stim Type=' num2str(ST) ] )
                    subplot(16,2,19)
                    plot(time,avgSe(10,:))
                    hold on
                    plot(time,avgDe(10,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch10, B=stand, R=Dev' ] )
                    subplot(16,2,20)
                    plot(time,avgSm(10,:))
                    hold on
                    plot(time,avgDm(10,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch10, B=stand, R=Dev' ] )
                      subplot(16,2,21)
                    plot(time,avgSe(11,:))
                    hold on
                    plot(time,avgDe(11,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch11, B=stand, R=Dev' ] )
                    subplot(16,2,22)
                    plot(time,avgSm(11,:))
                    hold on
                    plot(time,avgDm(11,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch11, B=stand, R=Dev' ] )
                      subplot(16,2,23)
                    plot(time,avgSe(12,:))
                    hold on
                    plot(time,avgDe(12,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch12, B=stand, R=Dev' ] )
                    subplot(16,2,24)
                    plot(time,avgSm(12,:))
                    hold on
                    plot(time,avgDm(12,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch12, B=stand, R=Dev' ] )
                      subplot(16,2,25)
                    plot(time,avgSe(13,:))
                    hold on
                    plot(time,avgDe(13,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch13, B=stand, R=Dev' ] )
                    subplot(16,2,26)
                    plot(time,avgSm(13,:))
                    hold on
                    plot(time,avgDm(13,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch13, B=stand, R=Dev' ] )
                    subplot(16,2,27)
                    plot(time,avgSe(14,:))
                    hold on
                    plot(time,avgDe(14,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch14, B=stand, R=Dev' ] )
                    subplot(16,2,28)
                    plot(time,avgSm(14,:))
                    hold on
                    plot(time,avgDm(14,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch14, B=stand, R=Dev' ] )
                    subplot(16,2,29)
                    plot(time,avgSe(15,:))
                    hold on
                    plot(time,avgDe(15,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch15, B=stand, R=Dev' ] )
                    subplot(16,2,30)
                    plot(time,avgSm(15,:))
                    hold on
                    plot(time,avgDm(15,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch15, B=stand, R=Dev' ] )
                    subplot(16,2,31)
                    plot(time,avgSe(16,:))
                    hold on
                    plot(time,avgDe(16,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'LFP Ch16. B=stand, R=Dev' ] )
                    subplot(16,2,32)
                    plot(time,avgSm(16,:))
                    hold on
                    plot(time,avgDm(16,:),'r')
                   set(gca,'xlim',[-25 200])
                    title( [ 'MUA Ch16, B=stand, R=Dev' ] )

                        axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,['2-' filenamesout  '-'  ', ' num2str(length(trig01)) ' sweeps'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                     print( '-djpeg','-r300',[directory1 '2-LFP__MUA_Stnd_Dev Stim Type=' num2str(ST) '_' filenamesout  '.jpg'])

avgDe2=avgDe;
avgDm2=avgDm;
avgSe2=avgSe;
avgSm2=avgSm;
XB=X;

 save([directory1 filenamesout], 'XA','XB','avgDe1','avgDm1','avgSe1','avgSm1','time','adrate','ss','D','trig01','avgDe2','avgDm2','avgSe2','avgSm2','newadrate','-mat','-v7.3')
end