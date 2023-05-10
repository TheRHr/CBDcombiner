%% Data combiner 
%All 3 loops are limited atm to only do the first item in directory
%to unlock, remove the '-5' or '-10' at the start of each for-loop
%requires Buzaki functions in "FeatureSynthesis" to be added to path
%requires "LFP_CBD_chronic_downsampled" to be in path directory
%requires a preformatted conditions table in path directory

clear all;

%set path from basepath to the dataset directory
path = "LFP_CBD_chronic_downsampled/";
directories = dir(path);

%read the pre-simplefied / formatted conditions table
CondTable = readtable('RatConditionsSimple.xlsx');


%initialize empty matrix to append to
X = [];

%set studydaycounter to 1 (for conditions table)
xCount = 0;

colnames = {'Ratnumber','StudyDay','posttrial (0 = presleep)','Condition (1=OR, 2=OD, 3=HC)','Treatment (0=VEH, 1=CBD)','Epoch (seconds)','DeltaHPC','DeltaPFC','ThetaHPC','ThetaPFC','EMG'};
samplingrate = 2500; 
TargetSampling = 1250;
timesDownSamp  = samplingrate / TargetSampling;
samplingFrequencyEMG = 5;
smoothWindowEMG = 10;

%loop over Rat-ID directories
for i =3:length(directories)-5
    ratNum = directories(i).name;
    %or just take i as ratNum...
    ratPath = append(path,directories(i).name,"/");
    ratDirs = dir(ratPath);

    %loop over studydays directories 
    for j =3:length(ratDirs)-10
        SD = ratDirs(j).name;
        SDPath = append(ratPath,ratDirs(j).name,"/");
        SDdirs = dir(SDPath);
        
        %get the studyday number from the name of the folder
        SDNum = str2double(extractBefore(extractAfter(SD,'SD'), '_'));
        
        %set counter to next studyday entry in the conditions table
        xCount = xCount + 1;
        
        %set trialcounter at 0 because first trial is presleep
        trialNum = 0;
        
        %loop over trial directories (presleep, posttrial1-5)
        for k = 3:length(SDdirs)-5
            trial = SDdirs(k).name;
            trialPath = append(SDPath,SDdirs(k).name,"/");
            trialDir = dir(trialPath);


            %load in the raw trial data
            HPCload = load(append(trialPath,"HPC_100_CH18_0.continuous.mat"));
            HPC =HPCload.HPC;
            PFCload = load(append(trialPath,"PFC_100_CH22_0.continuous.mat"));
            PFC =PFCload.PFC;

            %downsample the raw data from 2500 to 1250
            lfpPFCDown = decimate(PFC,timesDownSamp,'FIR');
            lfpHPCDown = decimate(HPC,timesDownSamp,'FIR');
            
            %Buzaki methods to create data for the feature matrix
            %identical to FeatureSynthesis
            timVect = linspace(0,numel(lfpPFCDown)/TargetSampling,numel(lfpPFCDown));
            DeltaBandPFC = compute_delta_buzsakiMethod(lfpPFCDown,timVect,TargetSampling,'DeltaBandPFCMat');
            DeltaBandHPC = compute_delta_buzsakiMethod(lfpHPCDown,timVect,TargetSampling,'DeltaBandHPCMat');
            ThetaBandPFC = compute_theta_buzsakiMethod(lfpPFCDown,timVect,TargetSampling,'ThetaBandPFCMat');
            ThetaBandHPC = compute_theta_buzsakiMethod(lfpHPCDown,timVect,TargetSampling,'ThetaBandHPCMat');
            BetaBandPFC = compute_beta_buzsakiMethod(lfpPFCDown,timVect,TargetSampling,'BetaBandPFCMat');
            BetaBandHPC = compute_beta_buzsakiMethod(lfpHPCDown,timVect,TargetSampling,'BetaBandHPCMat');
            GammaBandPFC = compute_gamma_buzsakiMethod(lfpPFCDown,timVect,TargetSampling,'GammaBandPFCMat');
            GammaBandHPC = compute_gamma_buzsakiMethod(lfpHPCDown,timVect,TargetSampling,'GammaBandHPCMat');
            EMGFromLFP = compute_emg_buzsakiMethod(samplingFrequencyEMG, TargetSampling, lfpPFCDown, lfpHPCDown, smoothWindowEMG,'EMGLikeSignalMat');
            prEMGtime = DeltaBandPFC.timestamps<EMGFromLFP.timestamps(1) | DeltaBandPFC.timestamps>EMGFromLFP.timestamps(end);
            DeltaBandPFC.data(prEMGtime) = []; 
            DeltaBandHPC.data(prEMGtime) = [];
            ThetaBandPFC.data(prEMGtime) = [];
            ThetaBandHPC.data(prEMGtime) = []; 
            GammaBandPFC.data(prEMGtime) = [];
            GammaBandHPC.data(prEMGtime) = [];
            DeltaBandPFC.timestamps(prEMGtime) = [];
            EMG = interp1(EMGFromLFP.timestamps,EMGFromLFP.smoothed,DeltaBandPFC.timestamps,'nearest');
            EMG = bz_NormToRange(EMG,[0 1]);

            %Initialize feature matrix and fill in all the columns
            lfpFeatures = zeros(length(EMG),5);
            lfpFeatures(:,1) = str2double(ratNum);
            lfpFeatures(:,2) = SDNum;
            lfpFeatures(:,3) = trialNum;  
            lfpFeatures(:,4) = table2array(CondTable(xCount,3));
            lfpFeatures(:,5) = table2array(CondTable(xCount,4));
            lfpFeatures(:,6) = 1:1:length(EMG);
            lfpFeatures(:,7) = DeltaBandHPC.data;
            lfpFeatures(:,8) = DeltaBandPFC.data;
            lfpFeatures(:,9) = ThetaBandHPC.data;
            lfpFeatures(:,10) = ThetaBandPFC.data;
            lfpFeatures(:,11) = EMG;

            %Append trial feature matrix to the global matrix 'X'
            X = [X ;lfpFeatures];
            
            %Trialcounter goes up for next trial
            trialNum = trialNum + 1;
        end

    end
  
end

%convert combined matrix to table format with the preset headers
X = array2table(X,'VariableNames',colnames);

%save combined table in both excel and mat file format
matfilename = 'CombinedCBDtable.mat';
save(matfilename,"X");
filename = 'CombinedCBDtable.xlsx';
writetable(X,filename);
