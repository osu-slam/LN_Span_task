%% Letter-Number Span Task v1
%
% The task involves the auditory presentation by an examiner of a mixed
% series of alternating numbers and letters (Gold et el., 1997).
% A subject is asked to respond by first saying the numbers in order from
% the smallest to the largest, followed by saying the letters in
% alphabetical order. For example, a subject would hear "w7t4," and
% the correct answer is "47tw."
% 
% To administer, press Run, fill out subject information, 
% and guide the participant through the instructions.
%
% Author -- Hyun-Woong Kim, The Ohio State University, khw173@gmail.com

sca; %DisableKeysForKbCheck([]); KbQueueStop;
clc; clear;

try 
    PsychPortAudio('Close'); 
catch
    disp('PsychPortAudio is already closed.')
end

testing = 0; % Change this to change expt window size.
textSize = 55; % Change this to change size of text. 


%% Collect subject information

rootdir = pwd;
expName = 'LN_span_task_v1';					    
expDate = date;

prompt = { 'Subject ID (only number!):', ...
           'Subject Information:', ...
           'Set A (1) or Set B (2):', ...
           'Skip Practice (1 to skip):'};
dlg_ans = inputdlg(prompt);
p.subjID  = str2double(dlg_ans{1});
p.subjInitial  = upper(dlg_ans{2});

if str2double(dlg_ans{3}) == 1
    p.session = 'pre';  p.stimSet = 'A';
elseif str2double(dlg_ans{3}) ==  2
    p.session = 'post';  p.stimSet = 'B';
else
    error('Please indicate if this is the pre- or post-testing session by entering 1 or 2')
end 

NoPractice = str2double(dlg_ans{4});
if NoPractice ~= 1
    NoPractice = 0;
end

% subject ID should be an integer number for between-subject conditions
if p.subjID~=round(p.subjID)
    error('invalid subject ID');
end

filename = strcat(rootdir, '/data/', expName, '_', num2digit(p.subjID,2), '_', p.session);
stimdir = [pwd '/stim_LN'];

DATAFOLDER = 'data';
if (~exist(DATAFOLDER,'dir'))
    mkdir(DATAFOLDER);
end

if exist([filename, '.mat'],'file')
    error('existing file name');
end


%% Set parameters

% setting variables using parameters above
p.sampleRate = 44100;  % sample rate
p.setSize = 2:7;  % set size list
p.trialPerSize = 4;  % number of trials per each set
p.intervalDuration = 1.5;  % interval duration between elements 

srate = p.sampleRate;  tps = p.trialPerSize;  intv=p.intervalDuration;
setSize = p.setSize;  nsize = length(p.setSize);

% setting variables using parameters above
L = ["C","D","F","H","J","K","L","P","Q","R","S","T","W"];
nL = length(L);  p.Letter = L;
N = strings(1,9);  N(1:9) = 1:9;
nN = length(N);  p.Number = N;

prac_A = ["D6","W4H"];  p.prac_A = prac_A;
prac_B = ["F8","R3K"];  p.prac_B = prac_B;

SP2_A = ["C3", "F2", "K8", "R7"];
SP3_A = ["D6T", "T9F", "P4R", "K5H"];
SP4_A = ["L3Q6","H4D7","F5J2","W7T4"];
SP5_A = ["4H6Q9","D8P5L","T3K6H","7J2F4"];
SP6_A = ["H8L4S1","K3F7R5","7Q3J2D","2W7D9J"];
SP7_A = ["1J3R5K8","D2P6J3T","F4W7Q9L","6C2S3F7"];
set_A = [SP2_A; SP3_A; SP4_A; SP5_A; SP6_A; SP7_A];
p.set_A = set_A;

SP2_B = ["T4",  "P9", "W5", "C6"];
SP3_B = ["S3W", "Q9D",  "L8P", "H6J"];
SP4_B = ["C4F7","K8P6","P5L3","J5H8"];
SP5_B = ["3P7R9","8K3H5","F7R5P","S4L6D"];
SP6_B = ["P3H6T4","1T5C8L","J7Q4W2","9R5K3F"];
SP7_B = ["2H4T6R9","D3S6J7H","7F1P4L8","C1L5F4W"];
set_B = [SP2_B; SP3_B; SP4_B; SP5_B; SP6_B; SP7_B];
p.set_B = set_B;


%% Stimulus conditions

p.condLabel = {'trialIndex','setSizeIndex','Letter-Number','correctness'};

if p.stimSet=='A'
    curr_prac = prac_A;  curr_set = set_A;
elseif p.stimSet=='B'
    curr_prac = prac_B;  curr_set = set_B;
end

% for practice block
p.numPracTrial = 2;
nprac = p.numPracTrial;
pracCond = strings(nprac,length(p.condLabel));
pracCond(:,1) = 1:p.numPracTrial;
pracCond(:,2) = 2:3;
pracCond(:,3) = curr_prac;

% for main expt
ntrial = nsize*tps;
cond = strings(ntrial,length(p.condLabel));
cond(:,1) = 1:ntrial;
for i=1:nsize
    idf = (i-1)*tps+1:i*tps;
    cond(idf,2) = setSize(i);
    cond(idf,3) = curr_set(i,:);
end


%% Generate or load auditory stimuli
% Same as preallocating variables, loading stimuli into Matlab before
% running the code helps keep Matlab's timing accurate. 

% generate audio outputs for Letters
audio_L = cell(nL,1);
for i=1:nL
    stimfile = [stimdir '/audio_' char(L(i)) '.wav'];
    [audio_tmp,fs] = audioread(stimfile);
    
    if fs~=srate  % double-check the sampling rate of audio files
        error('sampling rate is not equal to what you set');
    end
    audio_L{i} = audio_tmp';
end

% generate audio outputs for Numbers
audio_N = cell(nN,1);
for i=1:nN
    stimfile = [stimdir '/audio_' char(N(i)) '.wav'];
    [audio_tmp,fs] = audioread(stimfile);
    
    if fs~=srate  % double-check the sampling rate of audio files
        error('sampling rate is not equal to what you set');
    end
    audio_N{i} = audio_tmp';
end

% SPEAKER ICON AND FIXATION CROSS
speaker_mat = imread(fullfile(rootdir, 'Speaker_Icon.png'));
crossCoords = [-20, 20, 0, 0; 0, 0, -20, 20]; 


%% Open PsychToolbox (PTB) and RTBox
% PTB is used to generate the screen which the participant will see, and to
% present the auditory stimuli. If anyone is interested in learning to use 
% this incredibly powerful toolbox, I highly recommend checking out these 
% tutorials: http://peterscarfe.com/ptbtutorials.html
% [wPtr, rect] = Screen('OpenWindow', 0, 0);
Screen('Preference', 'SkipSyncTests', 1); 
if testing==1
    [w, rect] = Screen('OpenWindow', 0, 0, [0 0 800 600]);
else
    [w, rect] = Screen('OpenWindow', 0, 0);
end
p.ScreenIFI = Screen('GetFlipInterval', w);
cx = rect(3)/2;  cy = rect(4)/2;
Screen('TextSize', w, textSize);

DrawFormattedText(w, 'Please wait...', 'center', 'center', 255);
Screen('Flip', w);
WaitSecs(2);

ListenChar(2);
HideCursor(); 

InitializePsychSound(1);
pahandle = PsychPortAudio('Open', 1, [], 1, srate, 1);
% pahandle = PsychPortAudio('Open', 1, [], [], fs);

% RTBox is used to collect subject response and maintain timing of the
% experiment. It was originally designed for use in MRI, but I prefer to
% use it in behavioral experiments as well. There are very few tutorials
% online, so I recommend reading RTBox.m and RTBoxdemo.m 
RTBox('fake', 1);
RTBox('UntilTimeout', 1);
RTBox('ButtonNames', {'left', 'right', 'space', '4'});

% I convert the speaker image matrix into a texture at this point so the
% experiment runs faster. 
speaker_tex = Screen('MakeTexture', w, speaker_mat);


%% Practice block

if ~NoPractice
    % Instruction
    str = 'Welcome to the Letter-Number test.';
    str = [str '\n\nYou will hear a series of alternating numbers and letters.'];
    str = [str '\nPlease say the numbers in order from the smallest to the largest,'];
    str = [str '\nand then say the letters in alphabetical order.'];
    str = [str '\n\nFor example, if you hear "J3D1",'];
    str = [str '\nthe correct answer is "13DJ".'];
    str = [str '\n\nPress the space bar to begin a short practice'];

    DrawFormattedText(w, str, 'center', 'center', 255, [], [], [], 1.1);
    Screen('Flip', w);
    RTBox('Clear');
    RTBox(inf);
    
    Screen('Flip', w);
    WaitSecs(1);
    
    for t=1:nprac
        t_size = str2double(pracCond(t,2));
        t_str = char(pracCond(t,3));
        
        % auditory stimulus preparation
        audio_tmp = cell(t_size,1);
        for i=1:t_size
            iN = find(t_str(i)==N);  iL = find(t_str(i)==L);
            if ~isempty(iN)
                audio_tmp{i} = audio_N{iN};
            elseif ~isempty(iL)
                audio_tmp{i} = audio_L{iL};
            else
                error('string is incorrect');
            end
        end
        
        Screen('DrawLines', w, crossCoords, 2, 255, [cx, cy]);
        Screen('Flip', w);
        WaitSecs(1);
        
        Screen('DrawTexture', w, speaker_tex);
        Screen('Flip', w);
        for i=1:t_size
            PsychPortAudio('FillBuffer', pahandle, audio_tmp{i});
            PsychPortAudio('Start', pahandle);
            WaitSecs(intv);
        end
        Screen('Flip', w);
        RTBox('Clear'); 
        [~, answer] = RTBox(Inf);
        
        WaitSecs(1);
        if t<nprac
            str = 'Press the space bar to continue.';
            DrawFormattedText(w, str, 'center', 'center', 255, [], [], [], 1.1);
            Screen('Flip', w);
            WaitTill('space');
        end
    end
end


%% Actual Experiment

% Instruction
str = 'Great job! Let''s move on.';
str = [str '\n\nYou will hear a series of alternating numbers and letters.'];
str = [str '\nPlease say the numbers in order from the smallest to the largest,'];
str = [str '\nand then say the letters in alphabetical order.'];
str = [str '\n\nPress the space bar to begin.'];

DrawFormattedText(w, str, 'center', 'center', 255, [], [], [], 1.1);
Screen('Flip', w);
RTBox('Clear');
RTBox(inf);

Screen('Flip', w);
WaitSecs(1);

t=1;  quitit=0;  count=0;
while t<ntrial && ~quitit 
    t_size = str2double(cond(t,2));
    t_str = char(cond(t,3));

    % auditory stimulus preparation
    audio_tmp = cell(t_size,1);
    for i=1:t_size
        iN = find(t_str(i)==N);  iL = find(t_str(i)==L);
        if ~isempty(iN)
            audio_tmp{i} = audio_N{iN};
        elseif ~isempty(iL)
            audio_tmp{i} = audio_L{iL};
        else
            error('string is incorrect');
        end
    end

    Screen('DrawLines', w, crossCoords, 2, 255, [cx, cy]);
    Screen('Flip', w);
    WaitSecs(1);

    Screen('DrawTexture', w, speaker_tex);
    Screen('Flip', w);
    for i=1:t_size
        PsychPortAudio('FillBuffer', pahandle, audio_tmp{i});
        PsychPortAudio('Start', pahandle);
        WaitSecs(intv);
    end
    Screen('Flip', w);
    RTBox('Clear'); 
    [~, answer] = RTBox(Inf);

    if strcmp('left', answer),       cond(t,4)=0;
    elseif strcmp('right', answer),  cond(t,4)=1;  count=count+1;
    end

    WaitSecs(1);
    if t<ntrial
        str = 'Press the space bar to continue.';
        DrawFormattedText(w, str, 'center', 'center', 255, [], [], [], 1.1);
        Screen('Flip', w);
        key = WaitTill('space');
    end
    
    if mod(t,4)==0
        if count==0, quitit=1;  end
        count = 0;
    end
    
    % move onto next trial
    t=t+1;
end

if quitit==1
    p.highestStep = str2double(cond(t,2));
    idx = strcmp(cond(:,4),"");
    cond(idx,4) = 0;
end

% End of experiment
DrawFormattedText(w, 'End of the LN task.', 'center', 'center', 255);
Screen('Flip', w);
WaitSecs(4);

sca
ShowCursor; ListenChar();
save(filename, 'p','cond');