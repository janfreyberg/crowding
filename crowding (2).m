%% 
% made the targets slightly smaller (0.5 deg visual angle)
function crowding
try    
%% Get basic info, set filename.
rng('shuffle');
% subject info and screen info

ID = input('Participant ID? ', 's');
scr_diagonal = input('Screen Diagonal? ');
scr_distance = 57;
diagnosis = input('Diagnosis? ');

tstamp = clock;
if ~isdir( fullfile(pwd, 'Results', mfilename, num2str(diagnosis)) )
    mkdir( fullfile(pwd, 'Results', mfilename, num2str(diagnosis)) );
end
savefile = fullfile(pwd, 'Results', mfilename, num2str(diagnosis), [sprintf('crowding-%02d-%02d-%02d-%02d%02d-', tstamp(1), tstamp(2), tstamp(3), tstamp(4), tstamp(5)), ID, '.mat']);


%% Experiment Variables.
scr_background = 127.5;
scr_no = 1; %max(Screen('Screens'));
scr_dimensions = Screen('Rect', scr_no);
xcen = scr_dimensions(3)/2;
ycen = scr_dimensions(4)/2;

contrast = 1;
stimsize = 0.85;
stimdur = 2/60 - 0.01;
iti = 1;
isi = 3/60;

stimoffset = 10;
distances =  [0, 0.2, 0.4, 0.6, 1, 1.4];

questdist = 1;
questperf = 0.75;


% Set up Keyboard
WaitSecs(1);
KbName('UnifyKeyNames');
u_key = KbName('UpArrow');
d_key = KbName('DownArrow');
esc_key = KbName('Escape');
ent_key = KbName('Return'); ent_key = ent_key(1);
keyList = zeros(1, 256);
keyList([u_key, d_key, esc_key, ent_key]) = 1;
KbQueueCreate([], keyList); clear keyList


% Sound
InitializePsychSound;
pa = PsychPortAudio('Open', [], [], [], [], [], 256);
bp400 = PsychPortAudio('CreateBuffer', pa, [MakeBeep(400, 0.2); MakeBeep(400, 0.2)]);
PsychPortAudio('FillBuffer', pa, bp400);

% Open Window
scr = Screen('OpenWindow', scr_no, scr_background);
HideCursor;
if IsLinux
    Screen('TextFont', scr, '-schumacher-clean-medium-r-normal--0-0-75-75-c-0-koi8-r');
end
Screen('BlendFunction', scr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Priority(1);




%% Calculate Stimuli
stimoffset = visual_angle2pixel(stimoffset, scr_diagonal, scr_distance, scr_no);

% size in pix
stimsize = visual_angle2pixel(stimsize, scr_diagonal, scr_distance, scr_no);

fixpoint = [xcen, ycen];
fixlines = [-stimsize/3 stimsize/3  0 0;
             0 0 -stimsize/3 stimsize/3];
fixwidth = 6;

stimbox = [xcen-stimoffset-stimsize/2, xcen+stimoffset-stimsize/2
           ycen-stimsize/2,            ycen-stimsize/2
           xcen-stimoffset+stimsize/2, xcen+stimoffset+stimsize/2
           ycen+stimsize/2,            ycen+stimsize/2];

dbox = zeros(4, size(distances, 2), 2, 2);
i = 1;
for j = distances
    sep = stimsize + visual_angle2pixel(j, scr_diagonal, scr_distance, scr_no);
    dbox(:, i, 1, 1) = stimbox(:, 1) - [0; sep; 0; sep];
    dbox(:, i, 1, 2) = stimbox(:, 1) + [0; sep; 0; sep];
    dbox(:, i, 2, 1) = stimbox(:, 2) + [0; sep; 0; sep];
    dbox(:, i, 2, 2) = stimbox(:, 2) - [0; sep; 0; sep];
    i = i+1;
end

% gabor Patch
[x, y] = meshgrid(linspace(-3*pi, 3*pi, stimsize)); 
bell = (exp( -(x.^2/((3*pi/2)^2))-(y.^2/((3*pi/2)^2)) ));
gabor = (127.5*cos(y)).*contrast;
gabor = cat(3, 127.5+gabor, 127.5+gabor, 127.5+gabor, bell*255);
gabor = Screen('MakeTexture', scr, gabor);
clear x y bell



%% Calibrate the tobii scanner
% tetio_init;
% tetio_connectTracker('TT120-204-82700356');
% tetio_setFrameRate(120);
% PsychTobiiMirrorEyes(scr);
% PsychTobiiCalibrate(scr);




%% Demonstrate

Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
Screen('Flip', scr);
KbStrokeWait;

Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
Screen('DrawTexture', scr, gabor, [], stimbox(:, 1), 25);
Screen('Flip', scr);
KbStrokeWait;

Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
Screen('DrawTexture', scr, gabor, [], stimbox(:, 2), -25);
Screen('Flip', scr);
KbStrokeWait;

Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
Screen('DrawTexture', scr, gabor, [], stimbox(:, 1), 25);
Screen('DrawTexture', scr, gabor, [], dbox(:, 1, 1, 1));
Screen('DrawTexture', scr, gabor, [], dbox(:, 1, 1, 2));
Screen('Flip', scr);
KbStrokeWait;

Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
Screen('Flip', scr);
KbStrokeWait;




%% Practice 
% Hold
for ktrial = 1:20
    [practice.hold.correct(ktrial), practice.hold.rt(ktrial), practice.hold.gaze(ktrial)] = present(randi(2), 45*sign(rand - 0.5), 'hold', 1, randi(5));
    if practice.hold.correct(ktrial) < 0
        break
    end
end

% if mean(practice.hold.correct(11:20)) < 0.7 && sum(practice.hold.correct <0) == 0
%     error('Performance not high enough in the holding practice');
% end

break_img;

% Slow
[practice.slow.LR, practice.slow.tilt, practice.slow.distance] = make_schedule(1:6, 2);

practice.slow.n = size(practice.slow.LR, 2);

practice.slow.tilt = practice.slow.tilt .* (linspace(45, 20, practice.slow.n) + randi([-3, 3], 1, practice.slow.n));

practice.slow.dur = linspace(0.3, stimdur, practice.slow.n);

for ktrial = 1:practice.slow.n
    [practice.slow.correct(ktrial), practice.slow.rt(ktrial), practice.slow.gaze(ktrial)] = present(practice.slow.LR(ktrial), practice.slow.tilt(ktrial), practice.slow.dur(ktrial), 1, practice.slow.distance(ktrial));
    if practice.slow.correct(ktrial) < 0
        break
    end
end

if mean(practice.slow.correct) < 0.65 && sum(practice.slow.correct <0) == 0
    error('Performance not high enough in the slow practice');
end

break_img;

% Fast
[practice.fast.LR, practice.fast.tilt, practice.fast.distance] = make_schedule(1:6, 2);

practice.fast.n = size(practice.fast.LR, 2);

practice.fast.tilt = practice.fast.tilt .* (linspace(45, 20, practice.fast.n) + randi([-3, 3], 1, practice.fast.n));

for ktrial = 1:practice.fast.n
    [practice.fast.correct(ktrial), practice.fast.rt(ktrial), practice.fast.gaze(ktrial)] = present(practice.fast.LR(ktrial), practice.fast.tilt(ktrial), stimdur, 1, practice.fast.distance(ktrial));
    if practice.fast.correct(ktrial) < 0
        break
    end
end

if mean(practice.fast.correct) < 0.65 && sum(practice.fast.correct <0) == 0
    error('Performance not high enough in the fast practice');
end

break_img;




%% Staircases
order = randperm(2);
for stairi = 1:2
if order(stairi) == 1
%% Crowded Staircase

for i = 1:10
    [correct, ~] = present(randi(2), 35*sign(rand - 0.5), stimdur, 1, questdist);
    if correct < 0
        break
    end
end

[stair.crowded.LR, stair.crowded.tilt, stair.crowded.distance] = make_schedule(1, 13);

stair.crowded.tilt(1) = 35;
stair.crowded.reversals = 0;
n = size(stair.crowded.LR, 2);
for i = 1:n
    [stair.crowded.correct(i), stair.crowded.rt(i)] = present(stair.crowded.LR(i), stair.crowded.tilt(i), stimdur, 1, questdist);
    
    % break
    if stair.crowded.correct(i) < 0
        break
    end
    
    % reversals
    if i > 1 && stair.crowded.correct(i) ~= stair.crowded.correct(i-1)
        stair.crowded.reversals = stair.crowded.reversals + 1;
    end
    
    % update next tilt
    if i < size(stair.crowded.LR, 2)
        stair.crowded.tilt(i+1) = update_tilt(abs(stair.crowded.tilt(i)), stair.crowded.correct(i), i, stair.crowded.reversals) * sign(stair.crowded.tilt(i+1));
    end
    
end

stair.crowded.threshold = mean(abs(stair.crowded.tilt(n-9:n)));

break_img;

elseif order(stairi) == 2
%% Uncrowded staircase
for i = 1:10
    [correct, ~] = present(randi(2), 35*sign(rand - 0.5), stimdur, 0, questdist);
    if correct < 0
        break
    end
end

[stair.uncrowded.LR, stair.uncrowded.tilt, stair.uncrowded.distance] = make_schedule(1, 13);

stair.uncrowded.tilt(1) = 35;
stair.uncrowded.reversals = 0;
n = size(stair.uncrowded.LR, 2);
for i = 1:n
    [stair.uncrowded.correct(i), stair.uncrowded.rt(i)] = present(stair.uncrowded.LR(i), stair.uncrowded.tilt(i), stimdur, 0, questdist);
    
    % break
    if stair.uncrowded.correct(i) < 0
        break
    end
    
    % reversals
    if i > 1 && stair.uncrowded.correct(i) ~= stair.uncrowded.correct(i-1)
        stair.uncrowded.reversals = stair.uncrowded.reversals + 1;
    end
    
    % update next tilt
    if i < size(stair.uncrowded.LR, 2)
        stair.uncrowded.tilt(i+1) = update_tilt(abs(stair.uncrowded.tilt(i)), stair.uncrowded.correct(i), i, stair.uncrowded.reversals) * sign(stair.uncrowded.tilt(i+1));
    end
    
end

stair.uncrowded.threshold = mean(abs(stair.uncrowded.tilt(n-9:n)));

break_img;

end
end




%% Trial Blocks


[trials.LR, trials.tilt, trials.distance] = make_schedule(1:6, 9);
trials.tilt = trials.tilt *  stair.crowded.threshold;



for i = 1:size(trials.LR, 2)

    [trials.correct(i), trials.rt(i), trials.gaze(i)] = present(trials.LR(i), trials.tilt(i), stimdur, 1, trials.distance(i));
    
    if trials.correct(i) < 0
        break
    end
    
    % this is the stop check for performance at the THRESHOLDED distance, to see how well people are doing
    if sum(trials.distance(1:i) == 1) == 10 && abs(mean(trials.correct(trials.distance(1:i) == 1)) - 0.75) >= 0.1
        error('');
    end
    
    
    if i < size(trials.LR, 2) && mod(i, 54) == 0
        break_img;
    end
end




%% Shutdown
Screen('DrawText', scr, 'Thank you!', xcen-100, ycen);
Screen('Flip', scr);
WaitSecs(2);

save(savefile);
KbQueueStop;
PsychPortAudio('Close');
sca;
try
    tetio_cleanUp;
catch err2
    warning(err2.message)
end
Priority(0);




catch err
%% Catch
    KbQueueCheck;
    KbQueueStop;
    sca;
    PsychPortAudio('Close');
    savefile = [savefile(1:(size(savefile, 2)-4)), '-ERROR.mat'];
    save(savefile);
    try
        tetio_cleanUp;
    catch err2
        warning(err2.message)
    end
    if strcmp(input('Do you want to keep the data? y / n ', 's'), 'n')
        delete(savefile);
        disp('Data not saved.');
    end
    Priority(0);
    rethrow(err);

end

%% Subfunctions
function [correct, rt, gaze] = present(position, targettilt, targetduration, crowdedflag, distance)
    Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
    Screen('Flip', scr);
    while KbCheck;
    end
    if targettilt < 0
        correct_key = u_key;
        false_key = d_key;
    elseif targettilt > 0
        correct_key = d_key;
        false_key = u_key;
    elseif targettilt == 0
        if round(rand)
            correct_key = u_key;
            false_key = d_key;
        else
            correct_key = d_key;
            false_key = u_key;
        end
    end
    
    
    
    keyIsDown=0;
    WaitSecs(iti-0.5);
    tetio_startTracking;
    WaitSecs(0.5);
    
    
    
    if crowdedflag
        
        KbQueueStart;
        
        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        c_on = Screen('Flip', scr);

        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        Screen('DrawTexture', scr, gabor, [], stimbox(:, position), targettilt);
        Screen('DrawTexture', scr, gabor, [], dbox(:, distance, position, 1), 0);
        Screen('DrawTexture', scr, gabor, [], dbox(:, distance, position, 2), 0);
        t_on=Screen('Flip', scr, c_on+isi);
        
            if ischar(targetduration) && strcmp(targetduration, 'hold')
                targetduration = 0;
                while ~keyIsDown
                    [keyIsDown, first_press] = KbQueueCheck;
                    WaitSecs(0.05);
                end
            end

        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        t_off = Screen('Flip', scr, t_on+targetduration);

        while ~keyIsDown
            [keyIsDown, first_press] = KbQueueCheck;
            WaitSecs(0.05);
        end
        
% *************************************************************************
    else
        
        KbQueueStart;
        
        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        c_on = Screen('Flip', scr);

        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        Screen('DrawTexture', scr, gabor, [], stimbox(:, position), targettilt);
        t_on=Screen('Flip', scr, c_on+isi);
        
        if ischar(targetduration) && strcmp(targetduration, 'hold')
            while ~keyIsDown
                [keyIsDown, first_press] = KbQueueCheck;
            end
        end
                
        Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
        t_off = Screen('Flip', scr, t_on+targetduration);

        while ~keyIsDown
            [keyIsDown, first_press] = KbQueueCheck;
            WaitSecs(0.05);
        end
        
        
        
    end
    
    [gaze.lefteye, gaze.righteye, gaze.timestamp, ~] = tetio_readGazeData;
    tetio_stopTracking;
    
    
    if first_press(KbName('Escape'))
        error('You interrupted the script');
    elseif first_press(correct_key) && first_press(false_key) % if both keys are pressed, an error should be issued (this will virtually never happen)
        correct = 0;
        rt = 0;
        PsychPortAudio('Start', pa);
    elseif first_press(correct_key) % if correct, move on
        rt = first_press(correct_key)-t_on;
        correct = 1;
    elseif first_press(false_key)
        correct = 0;
        rt = first_press(false_key)-t_on;
        PsychPortAudio('Start', pa);
    elseif first_press(ent_key);
        rt = -1;
        correct = -1;
    end
    
    gaze.trial_onset = tetio_localToRemoteTime(int64(t_on * 10^6));
    gaze.response_moment = tetio_localToRemoteTime(int64((t_on + rt) * 10^6));
    
    KbQueueStop;
    if targetduration ~= 0 && abs(t_off - t_on) > targetduration+1/0.01
        warning('The timing in this trial was bad - presumably some flips were missed');
        disp({'stimdur', targetduration, 'actual dur', t_off - t_on});
    end
end

function break_img
    
    image = fullfile(pwd, 'Pictures', sprintf('%d.jpg', randi(50)));
    temp_image = imread(image, 'jpg');
    KbQueueStart;
    for h = 1:30
    msg = sprintf('Pause for %d s', 31-h);
    Screen('DrawText', scr, msg, xcen-100, 0);
    Screen('PutImage', scr, temp_image);
    Screen('Flip', scr);
    WaitSecs(1);
        [~, first_press] = KbQueueCheck;
        if first_press(esc_key)
            error('You interrupted the script by pressing Escape after exposure');
        elseif first_press(ent_key)
            break
        end
    end
    Screen('DrawLines', scr, fixlines, fixwidth, 255, fixpoint);
    Screen('Flip', scr);
    KbQueueStop;
    KbWait;
    WaitSecs(1);
end

function [new_tilt] = update_tilt(old_tilt, correct, ~, reversals)
if reversals <= 5
    stepsize = 0.15;
elseif reversals <= 10
    stepsize = 0.025;
else
    stepsize = 0.010;
end

if correct == 1
    new_tilt = old_tilt - old_tilt*stepsize *(1-questperf)*4;
elseif correct == 0
    new_tilt = old_tilt + old_tilt*stepsize *questperf*4;
end

if new_tilt > 80
    new_tilt = 80;
end

end

function [LR, tilt, distance] = make_schedule(distances, n)
    k = 0;
    for iLR = 1:2
        for iTilt = [-1, 1]
            for iDistance = distances
                for kk = 1:n
                k = k+1;
                distance(k) = iDistance;
                LR(k) = iLR;
                tilt(k) = iTilt;
                end
            end
        end
    end
    schedule = randperm(k);
    distance = distance(schedule);
    LR = LR(schedule);
    tilt = tilt(schedule);


end


end

