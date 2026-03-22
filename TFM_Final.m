clear all; close all; clc;
sep = '\';

% add ECGkit
addpath(genpath('C:\Users\cicir\Documents\TFM\ecg-kit-0.1.6'))

%% main folders
rootPa  = 'C:\Users\cicir\Documents\TFM\PatientSamples';
addpath(rootPa);

resultsDir = [rootPa sep 'Results after 2 years of treatment'];
ensure_dir(resultsDir);

%% ============================================================
% PART 1: PROCESS ALL ECGS, DETECT MARKERS, SAVE BIOMARKERS
%
% For every patient/date:
% 1) load BRLPF ECG
% 2) detect R peaks with PanTompkins
% 3) build median beat per lead
% 4) repeat median beat 10x -> synthetic signal
% 5) run delineation on synthetic signal
% 6) take beat 3
% 7) rescue missing P markers from tail of beat 2 if needed
% 8) compute biomarkers
% 9) save figures + biomarker mat file
%% ============================================================

ponLogRows = {}; % {Patient, DateFolder, Lead, Status, Marker_index, R_index, Note}
list_PA = dir(rootPa);

for idPa = 3:length(list_PA)

    if ~list_PA(idPa).isdir, continue; end
    NAME_pat = list_PA(idPa).name;
    if startsWith(NAME_pat,'.'), continue; end

    patDir = [rootPa sep NAME_pat];
    if ~exist(patDir,'dir'), continue; end

    list_Date = dir(patDir);

    for iDate = 3:length(list_Date)

        if ~list_Date(iDate).isdir, continue; end
        NAME_ecg = list_Date(iDate).name;
        if strcmp(NAME_ecg,'.') || strcmp(NAME_ecg,'..'), continue; end

        rec_folder = [patDir sep NAME_ecg];
        if ~exist(rec_folder,'dir'), continue; end

        %% load processed ECG + header
        signr = [NAME_pat '_' NAME_ecg 'BRLPF.mat'];
        sigFile = [rec_folder sep signr];
        if ~exist(sigFile,'file'), continue; end

        load(sigFile)                 % -> sinal
        signal = sinal';              % [N x 12]
        clear sinal

        hdrFile = [rec_folder sep NAME_pat '_' NAME_ecg '_HEADER.mat'];
        if ~exist(hdrFile,'file'), continue; end
        load(hdrFile)

        header = build_header_from_HEADER(HEADER, NAME_pat, NAME_ecg);
        labels = header.desc;
        if ~iscell(labels), labels = cellstr(labels); end
        labels = strtrim(labels);

        fs = header.freq;
        N  = size(signal,1);
        t_sec = (0:N-1)'/fs;

        baseName = erase_mat_ext(signr);

        %% save ECGkit-style copy + run PanTompkins on original signal
        SaveDIR_QRS = [rec_folder sep 'QRSDetection'];
        ensure_dir(SaveDIR_QRS);
        save([SaveDIR_QRS sep signr],'header','signal')

        ECG = ECGwrapper('recording_name',[SaveDIR_QRS sep signr]);
        ECG.ECGtaskHandle = 'QRS_detection';
        ECG.ECGtaskHandle.detectors = {'pantom'};
        ECG.Run;

        qrsFile = [SaveDIR_QRS sep NAME_pat '_' NAME_ecg 'BRLPF_QRS_detection.mat'];
        if ~exist(qrsFile,'file')
            clear ECG header HEADER
            continue;
        end
        load(qrsFile)

        %% plot all leads with PanTompkins R markers
        figAll = figure('Visible','off','Name',[NAME_pat ' ' NAME_ecg ' - all leads'],'NumberTitle','off');
        tiledlayout(4,3,"TileSpacing","compact","Padding","compact");

        for k = 1:12
            nexttile;
            plot(t_sec, signal(:,k)); hold on; grid on;
            title(labels{k});

            varname_p = ['pantom_' labels{k}];
            if exist(varname_p,'var')
                Pp = eval(varname_p);
                if isstruct(Pp) && isfield(Pp,'time') && ~isempty(Pp.time)
                    plot(t_sec(Pp.time), signal(Pp.time,k), 'k*', 'MarkerSize', 4);
                end
            end
            xlabel('s'); ylabel('mV');
        end
        save_figure_safe(figAll, [SaveDIR_QRS sep baseName '_AllLeads_Pantom.png']);

        %% build median beat per lead, window = R-250ms to R+750ms
        preR   = round(0.250 * fs);
        postR  = round(0.750 * fs);
        winLen = preR + postR + 1;
        tb     = (-preR:postR)'/fs;

        medBeats = cell(12,1);
        hasBeat  = false(12,1);

        for k = 1:12
            varname_p = ['pantom_' labels{k}];
            beatMat = [];

            if exist(varname_p,'var')
                Pp = eval(varname_p);
                if isstruct(Pp) && isfield(Pp,'time') && ~isempty(Pp.time)
                    rloc = Pp.time(:);
                    ok = rloc - preR >= 1 & rloc + postR <= N;
                    rloc = rloc(ok);

                    if ~isempty(rloc)
                        beatMat = nan(numel(rloc),winLen);
                        for b = 1:numel(rloc)
                            idx = (rloc(b)-preR):(rloc(b)+postR);
                            beatMat(b,:) = signal(idx,k).';
                        end
                    end
                end
            end

            if isempty(beatMat), continue; end

            medBeat = median(beatMat,1,'omitnan');
            medBeats{k} = medBeat;
            hasBeat(k) = true;

            figLead = figure('Visible','off','Name',sprintf('%s %s - %s beats', NAME_pat, NAME_ecg, labels{k}), 'NumberTitle','off');
            h1 = plot(tb, beatMat.', 'Color', [0.7 0.7 0.7 0.5]); hold on; grid on;
            h2 = plot(tb, medBeat, 'r', 'LineWidth', 2);
            xlabel('Time around R [s]'); ylabel('Amplitude [mV]');
            title(sprintf('%s - median heartbeat (pantom), n=%d', labels{k}, size(beatMat,1)));
            legend([h1(1) h2], {'Beats','Median'}, 'Location','best');

            save_figure_safe(figLead, [SaveDIR_QRS sep baseName '_' labels{k} '_MedianBeats_Pantom.png']);
        end

        if ~any(hasBeat)
            clear header HEADER ECG
            clearvars -regexp ^pantom_
            continue;
        end

        %% build synthetic 10-beat signal from median beats
        SaveDIR_Del  = [rec_folder sep 'Delineation']; ensure_dir(SaveDIR_Del);
        SaveDIR_Mark = [rec_folder sep 'Markers'];     ensure_dir(SaveDIR_Mark);

        numRep   = 10;
        nsampDel = winLen * numRep;
        signal_del = zeros(nsampDel,12);

        for k = 1:12
            if hasBeat(k)
                signal_del(:,k) = repmat(medBeats{k},1,numRep).';
            end
        end

        header_del = header;
        header_del.recname = [NAME_pat '_' NAME_ecg 'BRLPF_median10rep'];
        header_del.nsamp   = nsampDel;

        signrDel = [header_del.recname '.mat'];
        signal   = signal_del;
        header   = header_del; 
        save([SaveDIR_Del sep signrDel],'signal','header')

        %% delineation on synthetic signal
        ECG = ECGwrapper('recording_name',[SaveDIR_Del sep signrDel]);
        ECG.ECGtaskHandle = 'ECG_delineation';
        ECG.Run;

        baseNameDel = header_del.recname;
        resFileDel  = [SaveDIR_Del sep baseNameDel '_ECG_delineation.mat'];
        if ~exist(resFileDel,'file')
            clear header HEADER ECG header_del signal_del
            continue;
        end
        load(resFileDel); % wavedet
        t_del = (0:nsampDel-1)'/fs;

        %% PanTompkins on synthetic signal
        ECG_R = ECGwrapper('recording_name',[SaveDIR_Del sep signrDel]);
        ECG_R.ECGtaskHandle = 'QRS_detection';
        ECG_R.ECGtaskHandle.detectors = {'pantom'};
        ECG_R.Run;

        qrsSynFile = [SaveDIR_Del sep baseNameDel '_QRS_detection.mat'];
        if ~exist(qrsSynFile,'file')
            clear ECG_R
            continue;
        end
        load(qrsSynFile)

        %% initialize biomarker arrays
        QRSd_ms         = nan(1,12);
        PQ_ms           = nan(1,12);
        QT_ms           = nan(1,12);
        QRSampl_pp_mV   = nan(1,12);
        MaxPeak_mV      = nan(1,12);
        MinPeak_mV      = nan(1,12);
        Energ           = nan(1,12);
        Pavg            = nan(1,12);
        AreaQRS_abs     = nan(1,12);

        %% analyze beat 3 of synthetic signal
        beatNum = 3;
        beatIdx = (beatNum-1)*winLen + 1 : beatNum*winLen;

        for k = 1:12

            leadLabel = labels{k};

            % get lead delineation struct
            Wlead = get_wlead_struct(wavedet, leadLabel);

            % get synthetic R peaks for this lead
            varname_p_syn = ['pantom_' leadLabel];
            Rraw = [];
            if exist(varname_p_syn,'var')
                Ptmp = eval(varname_p_syn);
                if isstruct(Ptmp) && isfield(Ptmp,'time') && ~isempty(Ptmp.time)
                    Rraw = Ptmp.time(:);
                end
            end

            refr = round(0.200 * fs);
            Rclean = clean_r_list(Rraw, refr, nsampDel);

            expectedR = beatIdx(1) + preR;
            R3 = pick_r_for_beat(Rclean, expectedR, beatIdx, round(0.250*fs));

            if isempty(R3)
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, 'MISSING_R', NaN, NaN, 'No R found for beat3'};
                continue;
            end

            % get all markers for this lead
            [Pon_all, Ppk_all, Poff_all, Qon_all, Qoff_all, Ton_all, Tpk_all, Toff_all] = extract_markers(Wlead, nsampDel);

            % keep only beat 3 markers
            Pon3_c  = Pon_all( Pon_all  >= beatIdx(1) & Pon_all  <= beatIdx(end) );
            Ppk3_c  = Ppk_all( Ppk_all  >= beatIdx(1) & Ppk_all  <= beatIdx(end) );
            Poff3_c = Poff_all(Poff_all >= beatIdx(1) & Poff_all <= beatIdx(end) );

            Qon3_c  = Qon_all( Qon_all  >= beatIdx(1) & Qon_all  <= beatIdx(end) );
            Qoff3_c = Qoff_all(Qoff_all >= beatIdx(1) & Qoff_all <= beatIdx(end) );

            Ton3_c  = Ton_all( Ton_all  >= beatIdx(1) & Ton_all  <= beatIdx(end) );
            Tpk3_c  = Tpk_all( Tpk_all  >= beatIdx(1) & Tpk_all  <= beatIdx(end) );
            Toff3_c = Toff_all(Toff_all >= beatIdx(1) & Toff_all <= beatIdx(end) );

            % main markers around beat 3
            Qon  = pick_before(Qon3_c,  R3, round(0.20*fs), 0);
            Qoff = pick_after (Qoff3_c, R3, round(0.25*fs), 0);

            Pon  = pick_before(Pon3_c,  R3, round(0.40*fs), round(0.02*fs));
            Ppk  = pick_before(Ppk3_c,  R3, round(0.35*fs), round(0.01*fs));
            Poff = pick_before(Poff3_c, R3, round(0.30*fs), 0);

            %% rescue missing P markers from tail of beat 2
            ponStatus  = 'OK';
            ppkStatus  = 'OK';
            poffStatus = 'OK';

            beat2Idx  = (beatNum-2)*winLen + 1 : (beatNum-1)*winLen;
            tailLen   = round(0.150*fs);
            tailStart = max(beat2Idx(end) - tailLen, beat2Idx(1));

            Pon2_tail  = Pon_all( Pon_all  >= tailStart & Pon_all  <= beat2Idx(end) );
            Ppk2_tail  = Ppk_all( Ppk_all  >= tailStart & Ppk_all  <= beat2Idx(end) );
            Poff2_tail = Poff_all(Poff_all >= tailStart & Poff_all <= beat2Idx(end) );

            prPon_min  = round(0.080*fs); prPon_max  = round(0.400*fs);
            prPpk_min  = round(0.060*fs); prPpk_max  = round(0.350*fs);
            prPoff_min = round(0.040*fs); prPoff_max = round(0.300*fs);

            Pon2_tail  = Pon2_tail(  Pon2_tail  <= (R3 - prPon_min)  & Pon2_tail  >= (R3 - prPon_max)  );
            Ppk2_tail  = Ppk2_tail(  Ppk2_tail  <= (R3 - prPpk_min)  & Ppk2_tail  >= (R3 - prPpk_max)  );
            Poff2_tail = Poff2_tail( Poff2_tail <= (R3 - prPoff_min) & Poff2_tail >= (R3 - prPoff_max) );

            if isnan(Pon)
                if ~isempty(Pon2_tail), Pon = Pon2_tail(end); ponStatus = 'RESCUED_FROM_BEAT2';
                else, ponStatus = 'MISSING'; end
            end
            if isnan(Ppk)
                if ~isempty(Ppk2_tail), Ppk = Ppk2_tail(end); ppkStatus = 'RESCUED_FROM_BEAT2';
                else, ppkStatus = 'MISSING'; end
            end
            if isnan(Poff)
                if ~isempty(Poff2_tail), Poff = Poff2_tail(end); poffStatus = 'RESCUED_FROM_BEAT2';
                else, poffStatus = 'MISSING'; end
            end

            % ordering checks for P markers
            if ~isnan(Pon) && ~isnan(Ppk)  && ~(Pon < Ppk),   Ppk  = NaN; ppkStatus  = 'DROPPED_BAD_ORDER'; end
            if ~isnan(Ppk) && ~isnan(Poff) && ~(Ppk < Poff),  Poff = NaN; poffStatus = 'DROPPED_BAD_ORDER'; end
            if ~isnan(Pon) && ~isnan(Poff) && ~(Pon < Poff),  Poff = NaN; poffStatus = 'DROPPED_BAD_ORDER'; end

            if ~isnan(Qon)
                if ~isnan(Poff) && ~(Poff < Qon), Poff = NaN; poffStatus = 'DROPPED_AFTER_QRS'; end
                if ~isnan(Ppk)  && ~(Ppk  < Qon), Ppk  = NaN; ppkStatus  = 'DROPPED_AFTER_QRS'; end
                if ~isnan(Pon)  && ~(Pon  < Qon), Pon  = NaN; ponStatus  = 'DROPPED_AFTER_QRS'; end
            end

            % log what happened with P markers
            if isnan(Pon)
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['PON_'  ponStatus], NaN, R3, 'No P_on found (beat3 nor tail beat2)'};
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['PPK_'  ppkStatus], Ppk, R3, ''};
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['POFF_' poffStatus], Poff, R3, ''};
                continue;
            else
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['PON_'  ponStatus], Pon,  R3, ''};
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['PPK_'  ppkStatus], Ppk,  R3, ''};
                ponLogRows(end+1,:) = {NAME_pat, NAME_ecg, leadLabel, ['POFF_' poffStatus], Poff, R3, ''};
            end

            %% pick T markers after QRS
            refT = R3;
            if ~isnan(Qoff), refT = Qoff; end

            Ton  = pick_after(Ton3_c,  refT, round(0.90*fs), round(0.02*fs));
            Tpk  = pick_after(Tpk3_c,  refT, round(1.10*fs), round(0.04*fs));
            Toff = pick_after(Toff3_c, refT, round(1.30*fs), round(0.08*fs));

            if ~(isnan(Ton) || isnan(Tpk))  && ~(Ton < Tpk),  Ton = NaN; Tpk = NaN; end
            if ~(isnan(Tpk) || isnan(Toff)) && ~(Tpk < Toff), Tpk = NaN; Toff = NaN; end

            %% plot beat 3 markers
            needTail = strcmp(ponStatus,'RESCUED_FROM_BEAT2') || strcmp(ppkStatus,'RESCUED_FROM_BEAT2') || strcmp(poffStatus,'RESCUED_FROM_BEAT2');
            if needTail, plotIdx = [tailStart:beat2Idx(end), beatIdx];
            else,        plotIdx = beatIdx; end
            plotIdx = unique(plotIdx);

            sigPlot = signal_del(plotIdx, k);
            tPlot   = t_del(plotIdx);
            t0      = t_del(Pon);
            t_rel   = tPlot - t0;

            figBeat = figure('Visible','off','NumberTitle','off');
            plot(t_rel, sigPlot); hold on; grid on;
            xlabel('Time relative to P_{on} [s]');
            ylabel('Amplitude [mV]');
            if needTail
                title(sprintf('%s - beat3 (P rescued: tail + beat3)', leadLabel), 'Interpreter','none');
            else
                title(sprintf('%s - 3rd beat (P_{on} anchored)', leadLabel), 'Interpreter','none');
            end

            leg2 = {'ECG'};
            plot(t_del(Pon)-t0,  signal_del(Pon,k),  'go', 'MarkerSize', 4); leg2{end+1}='P_{on}';
            if ~isnan(Ppk),  plot(t_del(Ppk)-t0,  signal_del(Ppk,k),  'g^', 'MarkerSize', 4); leg2{end+1}='P_{peak}'; end
            if ~isnan(Poff), plot(t_del(Poff)-t0, signal_del(Poff,k), 'gv', 'MarkerSize', 4); leg2{end+1}='P_{off}'; end
            if ~isnan(Qon),  plot(t_del(Qon)-t0,  signal_del(Qon,k),  'bo', 'MarkerSize', 4); leg2{end+1}='QRS_{on}'; end
            if ~isnan(Qoff), plot(t_del(Qoff)-t0, signal_del(Qoff,k), 'bx', 'MarkerSize', 4); leg2{end+1}='QRS_{off}'; end
            if ~isnan(Ton),  plot(t_del(Ton)-t0,  signal_del(Ton,k),  'mo', 'MarkerSize', 4); leg2{end+1}='T_{on}'; end
            if ~isnan(Tpk),  plot(t_del(Tpk)-t0,  signal_del(Tpk,k),  'm^', 'MarkerSize', 4); leg2{end+1}='T_{peak}'; end
            if ~isnan(Toff), plot(t_del(Toff)-t0, signal_del(Toff,k), 'mv', 'MarkerSize', 4); leg2{end+1}='T_{off}'; end
            plot(t_del(R3)-t0, signal_del(R3,k), 'r*', 'MarkerSize', 5); leg2{end+1}='R_{pantom}';
            legend(leg2,'Location','bestoutside');

            save_figure_safe(figBeat, [SaveDIR_Mark sep baseNameDel '_' leadLabel '_Beat3Markers.png']);

            %% compute biomarkers
            if ~(~isnan(Pon) && ~isnan(Qon) && ~isnan(Qoff) && ~isnan(Ton) && ~isnan(Tpk) && ~isnan(Toff))
                continue;
            end
            if ~(Qoff > Qon && Qon > Pon && Toff > Qon)
                continue;
            end

            QRSd_ms(1,k) = (Qoff - Qon) / fs * 1000;
            PQ_ms(1,k)   = (Qon  - Pon) / fs * 1000;
            QT_ms(1,k)   = (Toff - Qon) / fs * 1000;

            % local baseline before QRS
            n30     = round(0.030 * fs);
            startBL = max(Qon - n30, beatIdx(1));
            segPR   = median(signal_del(startBL:Qon,k));

            auxQRS = signal_del(Qon:Qoff,k) - segPR;
            if isempty(auxQRS), continue; end

            MaxPeak_mV(1,k)    = max(auxQRS);
            MinPeak_mV(1,k)    = min(auxQRS);
            QRSampl_pp_mV(1,k) = MaxPeak_mV(1,k) - MinPeak_mV(1,k);

            Energ(1,k)       = (1/fs) * sum(auxQRS.^2);
            Pavg(1,k)        = (1/numel(auxQRS)) * sum(auxQRS.^2);
            AreaQRS_abs(1,k) = (1/fs) * sum(abs(auxQRS)) * 1000;
        end

        %% save biomarkers to Markers folder
        biomarkers.leads          = labels;
        biomarkers.fs             = fs;
        biomarkers.QRSd_ms        = QRSd_ms;
        biomarkers.PQ_ms          = PQ_ms;
        biomarkers.QT_ms          = QT_ms;
        biomarkers.QRSampl_pp_mV  = QRSampl_pp_mV;
        biomarkers.MaxPeak_mV     = MaxPeak_mV;
        biomarkers.MinPeak_mV     = MinPeak_mV;
        biomarkers.Energ          = Energ;
        biomarkers.Pavg           = Pavg;
        biomarkers.AreaQRS_abs    = AreaQRS_abs;

        biomFile = [SaveDIR_Mark sep baseNameDel '_Beat3_Biomarkers.mat'];
        save(biomFile, 'biomarkers');

        %% cleanup per date
        clear header HEADER ECG ECG_R header_del signal_del wavedet biomarkers
        clearvars -regexp ^pantom_
    end
end

%% ================================
% save P-marker rescue log
%% ================================

PonT = cell2table(ponLogRows, 'VariableNames', {'Patient','DateFolder','Lead','Status','Marker_index','R_index','Note'});
ponXLSX = safe_writetable(PonT, [rootPa sep 'P_markers_Missing_or_Rescued.xlsx']);

%% =======================================================
% PART 2: BUILD COHORTS
%
% Rules:
% - Year1 = first usable folder in baseline year
% - Year3 = first usable folder in baseline year + 2
% - Year5 = first usable folder in baseline year + 4
% - usable means core markers exist in at least 8 leads
%% =======================================================

% QT excluded from Part 3 analysis
biomarkerFields_part3 = { ...
    'QRSd_ms', 'PQ_ms', 'QRSampl_pp_mV', ...
    'MaxPeak_mV','MinPeak_mV','Energ','Pavg','AreaQRS_abs'};

coreFields = {'QRSd_ms','PQ_ms','QT_ms'};
minValidLeads = 8;

selectionRows_13 = {};
selectionRows_15 = {};

includedPatients_13 = {};
includedY1_13 = {};
includedY3_13 = {};

includedPatients_15 = {};
includedY1_15 = {};
includedY5_15 = {};

list_PA = dir(rootPa);

for idPa = 3:length(list_PA)

    if ~list_PA(idPa).isdir, continue; end
    NAME_pat = list_PA(idPa).name;
    if startsWith(NAME_pat,'.'), continue; end

    patDir = [rootPa sep NAME_pat];
    if ~exist(patDir,'dir'), continue; end

    [dateNames, dateDT] = get_sorted_datefolders(patDir);

    if isempty(dateDT)
        selectionRows_13(end+1,:) = {NAME_pat, '', '', 'NoDateFolders'}; 
        selectionRows_15(end+1,:) = {NAME_pat, '', '', 'NoDateFolders', 'No'}; 
        continue;
    end

    baselineYear = year(dateDT(1));
    y1_year = baselineYear;
    y3_year = baselineYear + 2;
    y5_year = baselineYear + 4;

    % --- Year1 vs Year3
    [row13, included13, y1_date_13, y3_date] = select_single_comparison( ...
        rootPa, NAME_pat, dateNames, dateDT, y1_year, y3_year, ...
        'Year3', coreFields, minValidLeads, sep);

    selectionRows_13(end+1,:) = row13; 
    if included13
        includedPatients_13{end+1} = NAME_pat;
        includedY1_13{end+1} = y1_date_13; 
        includedY3_13{end+1} = y3_date; 
    end

    % --- Year1 vs Year5
    overlapTxt = ternary_str(ismember(NAME_pat, includedPatients_13), 'Yes', 'No');
    [row15base, included15, y1_date_15, y5_date] = select_single_comparison( ...
        rootPa, NAME_pat, dateNames, dateDT, y1_year, y5_year, ...
        'Year5', coreFields, minValidLeads, sep);

    row15 = [row15base, {overlapTxt}];
    selectionRows_15(end+1,:) = row15; 

    if included15
        includedPatients_15{end+1} = NAME_pat; 
        includedY1_15{end+1} = y1_date_15; 
        includedY5_15{end+1} = y5_date; 
    end
end

% refresh Y1vs5 overlap once Y1vs3 cohort is fully known
for i = 1:size(selectionRows_15,1)
    patName = selectionRows_15{i,1};
    selectionRows_15{i,5} = ternary_str(ismember(patName, includedPatients_13), 'Yes', 'No');
end

SelT13 = cell2table(selectionRows_13, 'VariableNames', {'Patient','Year1_DateFolder','Year3_DateFolder','Status'});
selectionLog13XLSX = safe_writetable(SelT13, [rootPa sep 'Year1_vs_Year3_SelectedFolders.xlsx']);

SelT15 = cell2table(selectionRows_15, 'VariableNames', {'Patient','Year1_DateFolder','Year5_DateFolder','Status','AlsoInYear1vs3'});
selectionLog15XLSX = safe_writetable(SelT15, [rootPa sep 'Year1_vs_Year5_SelectedFolders.xlsx']);

%% ============================================================
% PART 3: LONGITUDINAL ANALYSIS
%
% Use one reusable function for both:
% - Year1 vs Year3
% - Year1 vs Year5
%
% Global stats rule:
% if ANY biomarker is non-normal or lillietest fails,
% use Wilcoxon signrank for ALL biomarkers in that comparison
%% ============================================================

alphaLevel = 0.05;

[outTests13XLSX, forceWilcoxon13] = run_longitudinal_comparison( ...
    rootPa, resultsDir, sep, ...
    includedPatients_13, includedY1_13, includedY3_13, ...
    biomarkerFields_part3, ...
    'Year1_vs_Year3', 'Y1_vs_Y3', 'Year3', ...
    alphaLevel);

[outTests15XLSX, forceWilcoxon15] = run_longitudinal_comparison( ...
    rootPa, resultsDir, sep, ...
    includedPatients_15, includedY1_15, includedY5_15, ...
    biomarkerFields_part3, ...
    'Year1_vs_Year5', 'Y1_vs_Y5', 'Year5', ...
    alphaLevel, includedPatients_13);

%% final printout
fprintf('\nDone.\n');
fprintf('P-marker log: %s\n', ponXLSX);
fprintf('Selection log Y1vsY3: %s\n', selectionLog13XLSX);
fprintf('Selection log Y1vsY5: %s\n', selectionLog15XLSX);
fprintf('Paired tests Y1vsY3: %s\n', outTests13XLSX);
fprintf('Paired tests Y1vsY5: %s\n', outTests15XLSX);
fprintf('Per-patient plots saved under: %s\n', resultsDir);
fprintf('Global Wilcoxon forced for Y1vsY3: %s\n', ternary_str(forceWilcoxon13,'Yes','No'));
fprintf('Global Wilcoxon forced for Y1vsY5: %s\n\n', ternary_str(forceWilcoxon15,'Yes','No'));

%% =====================
% HELPERS
%% =====================

function header = build_header_from_HEADER(HEADER, NAME_pat, NAME_ecg)
    % rebuild the header in the format used later by ECGkit
    header.recname = [NAME_pat '_' NAME_ecg 'BRLPF'];
    header.desc    = HEADER.desc;
    header.freq    = HEADER.freq;
    header.nsamp   = HEADER.nsamp;
    header.nsig    = HEADER.nsig;
    header.gain    = HEADER.gain;
    header.adczero = HEADER.adczero;
    header.units   = HEADER.units;
    header.btime   = HEADER.btime;
    header.bdate   = HEADER.bdate;
end

function out = erase_mat_ext(fname)
    % remove .mat if it is there
    out = fname;
    if length(out) > 4 && strcmpi(out(end-3:end),'.mat')
        out = out(1:end-4);
    end
end

function ensure_dir(folderPath)
    % make directory if it does not exist
    if ~exist(folderPath,'dir'), mkdir(folderPath); end
end

function save_figure_safe(figHandle, outFile)
    % save figure using exportgraphics if possible, otherwise saveas
    try
        exportgraphics(figHandle, outFile);
    catch
        saveas(figHandle, outFile);
    end
    close(figHandle);
end

function outFile = safe_writetable(T, outFile)
    % write Excel table, if file is locked then save timestamped copy
    try
        writetable(T, outFile);
    catch
        [folderPath, baseName, ext] = fileparts(outFile);
        ts = datestr(now,'yyyymmdd_HHMMSS');
        outFile = [folderPath filesep baseName '_' ts ext];
        writetable(T, outFile);
    end
end

function safe_write_textfile(filePath, linesCell)
    % write text file from a cell array of lines
    try
        fid = fopen(filePath,'w');
        for i = 1:numel(linesCell)
            fprintf(fid,'%s\n', linesCell{i});
        end
        fclose(fid);
    catch
    end
end

function Wlead = get_wlead_struct(wavedet, leadLabel)
    % fetch lead struct from wavedet, handling upper/lower case too
    Wlead = [];
    if exist('wavedet','var') && isstruct(wavedet)
        if isfield(wavedet, leadLabel)
            Wlead = wavedet.(leadLabel);
        elseif isfield(wavedet, upper(leadLabel))
            Wlead = wavedet.(upper(leadLabel));
        elseif isfield(wavedet, lower(leadLabel))
            Wlead = wavedet.(lower(leadLabel));
        end
    end
end

function [dateNames, dateDT] = get_sorted_datefolders(patDir)
    % read patient date folders that look like yyyymmdd and sort them
    dList = dir(patDir);
    dateNames = {};
    dateDT = datetime.empty;

    for i = 1:length(dList)
        if ~dList(i).isdir, continue; end
        dn = dList(i).name;

        if length(dn)==8 && all(isstrprop(dn,'digit'))
            try
                dt = datetime(dn,'InputFormat','yyyyMMdd');
                dateNames{end+1} = dn; 
                dateDT(end+1) = dt; 
            catch
            end
        end
    end

    if ~isempty(dateDT)
        [dateDT, idxSort] = sort(dateDT);
        dateNames = dateNames(idxSort);
    end
end

function [rowOut, includedFlag, y1_date, yFollow_date] = select_single_comparison( ...
    rootPa, NAME_pat, dateNames, dateDT, y1_year, yFollow_year, ...
    followLabel, coreFields, minValidLeads, sep)
    % helper for Part 2:
    % builds one Year1-vs-followup group entry using the exact same logic

    includedFlag = false;
    y1_date = '';
    yFollow_date = '';

    y1_date = find_first_usable_folder_in_year(rootPa, NAME_pat, dateNames, dateDT, y1_year, coreFields, minValidLeads, sep);

    if isempty(y1_date)
        rowOut = {NAME_pat, '', '', 'NoUsableYear1'};
        return;
    end

    yFollow_date = find_first_usable_folder_in_year(rootPa, NAME_pat, dateNames, dateDT, yFollow_year, coreFields, minValidLeads, sep);

    if isempty(yFollow_date)
        rowOut = {NAME_pat, y1_date, '', ['NoUsable' followLabel]};
        return;
    end

    fileY1 = find_first_marker_biom_file(rootPa, NAME_pat, y1_date, sep);
    fileY2 = find_first_marker_biom_file(rootPa, NAME_pat, yFollow_date, sep);

    if isempty(fileY1) || isempty(fileY2)
        rowOut = {NAME_pat, y1_date, yFollow_date, 'MissingBiomarkerFile'};
        return;
    end

    rowOut = {NAME_pat, y1_date, yFollow_date, 'Included'};
    includedFlag = true;
end

function [outTestsXLSX, forceWilcoxon] = run_longitudinal_comparison( ...
    rootPa, resultsDir, sep, ...
    includedPatients, includedY1, includedY2, ...
    biomarkerFields_part3, ...
    comparisonLabel, suffixLabel, followupLabel, ...
    alphaLevel, varargin)
    % reusable longitudinal analysis block
    %
    % does:
    % 1) load Year1 + follow-up biomarker files
    % 2) save patient-specific plots
    % 3) collect patient-level medians
    % 4) run global stats table

    [Y1_med, Y2_med, D_med] = init_stat_structs(biomarkerFields_part3);

    overlapPatients = {};
    if ~isempty(varargin)
        overlapPatients = varargin{1};
    end

    for iP = 1:numel(includedPatients)

        pat = includedPatients{iP};
        y1d = includedY1{iP};
        y2d = includedY2{iP};

        fileY1 = find_first_marker_biom_file(rootPa, pat, y1d, sep);
        fileY2 = find_first_marker_biom_file(rootPa, pat, y2d, sep);
        if isempty(fileY1) || isempty(fileY2), continue; end

        S1 = load(fileY1); S2 = load(fileY2);
        if ~isfield(S1,'biomarkers') || ~isfield(S2,'biomarkers'), continue; end

        B1 = S1.biomarkers;
        B2 = S2.biomarkers;

        % patient result folder for this comparison
        patResDir = [resultsDir sep pat sep comparisonLabel];
        ensure_dir(patResDir);

        % save small text file showing which folders were used
        txtLines = {
            ['Patient: ' pat]
            ['Year1 folder: ' y1d]
            [followupLabel ' folder: ' y2d]
            };

        if ~isempty(overlapPatients)
            txtLines{end+1} = ['Also in Year1_vs_Year3 cohort: ' ternary_str(ismember(pat, overlapPatients), 'Yes', 'No')];
        end

        safe_write_textfile([patResDir sep pat '_' suffixLabel '_SelectedFolders.txt'], txtLines);

        % lead labels
        if isfield(B1,'leads')
            leadLabels = B1.leads;
            if ~iscell(leadLabels), leadLabels = cellstr(leadLabels); end
            leadLabels = strtrim(leadLabels);
        else
            leadLabels = {'I','II','III','aVL','aVR','aVF','V1','V2','V3','V4','V5','V6'};
        end

        % loop through biomarkers included in Part 3
        for f = 1:numel(biomarkerFields_part3)

            field = biomarkerFields_part3{f};
            if ~isfield(B1,field) || ~isfield(B2,field), continue; end

            v1 = B1.(field);
            v2 = B2.(field);
            if isempty(v1) || isempty(v2), continue; end

            v1 = v1(:).';
            v2 = v2(:).';
            if numel(v1) ~= 12 || numel(v2) ~= 12, continue; end

            % lead-wise difference = follow-up - Year1
            dLead = v2 - v1;

            % difference bar plot
            fig = figure('Visible','off','NumberTitle','off');
            bar(dLead(:));
            grid on;
            xticks(1:12);
            xticklabels(leadLabels);
            title(sprintf('%s | %s vs %s | %s (%s - Year1)', pat, y1d, y2d, field, followupLabel), 'Interpreter','none');
            ylabel_for_field(field);
            save_figure_safe(fig, [patResDir sep pat '_' y1d '_vs_' y2d '_' field '_DiffBars.png']);

            % boxplot of lead-wise differences
            fig2 = figure('Visible','off','NumberTitle','off');
            boxplot(dLead(:), 'Labels', {'Leads'});
            grid on;
            title(sprintf('%s | %s | Lead diffs (%s-Year1): %s', pat, [y1d ' -> ' y2d], followupLabel, field), 'Interpreter','none');
            ylabel_for_field(field);
            save_figure_safe(fig2, [patResDir sep pat '_' y1d '_vs_' y2d '_' field '_DiffBox.png']);

            % patient-level medians for stats
            m1 = median(v1,'omitnan');
            m2 = median(v2,'omitnan');

            if ~(isnan(m1) || isnan(m2))
                Y1_med.(field) = [Y1_med.(field); m1]; %
                Y2_med.(field) = [Y2_med.(field); m2]; 
                D_med.(field)  = [D_med.(field); (m2-m1)]; 
            end
        end
    end

    % build global test table
    [TestT, forceWilcoxon] = build_global_test_table( ...
        biomarkerFields_part3, Y1_med, Y2_med, D_med, alphaLevel);

    outTestsXLSX = safe_writetable(TestT, [rootPa sep comparisonLabel '_ParametricCheck_and_PairedTests.xlsx']);
end

function Rclean = clean_r_list(Rraw, refr, nsamp)
    % remove repeated / invalid R peaks and enforce refractory period
    Rclean = [];
    if isempty(Rraw), return; end

    Rraw = unique(Rraw(:));
    Rraw = Rraw(~isnan(Rraw) & Rraw >= 1 & Rraw <= nsamp);
    if isempty(Rraw), return; end

    Rclean = Rraw(1);
    last = Rraw(1);

    for i = 2:numel(Rraw)
        if Rraw(i) - last > refr
            Rclean(end+1,1) = Rraw(i); %#ok<AGROW>
            last = Rraw(i);
        end
    end
end

function R3 = pick_r_for_beat(Rclean, expectedR, beatIdx, margin)
    % pick the R that best matches the beat we want
    R3 = [];
    if isempty(Rclean), return; end

    lo = beatIdx(1) - margin;
    hi = beatIdx(end) + margin;
    cand = Rclean(Rclean >= lo & Rclean <= hi);

    if isempty(cand)
        [~,ii] = min(abs(Rclean - expectedR));
        R3 = Rclean(ii);
        return;
    end

    [~,ii] = min(abs(cand - expectedR));
    R3 = cand(ii);
end

function [Pon, Ppk, Poff, Qon, Qoff, Ton, Tpk, Toff] = extract_markers(Wlead, nsamp)
    % extract all possible marker arrays and keep only valid indices
    Pon=[]; Ppk=[]; Poff=[]; Qon=[]; Qoff=[]; Ton=[]; Tpk=[]; Toff=[];
    if isempty(Wlead), return; end

    clean = @(x) unique(x(~isnan(x) & x>=1 & x<=nsamp));

    if isfield(Wlead,'P_on'), Pon = clean(Wlead.P_on(:));
    elseif isfield(Wlead,'Pon'), Pon = clean(Wlead.Pon(:)); end

    if isfield(Wlead,'P'), Ppk = clean(Wlead.P(:));
    elseif isfield(Wlead,'Ppos'), Ppk = clean(Wlead.Ppos(:)); end

    if isfield(Wlead,'P_off'), Poff = clean(Wlead.P_off(:));
    elseif isfield(Wlead,'Poff'), Poff = clean(Wlead.Poff(:)); end

    if isfield(Wlead,'QRS_on'), Qon = clean(Wlead.QRS_on(:));
    elseif isfield(Wlead,'QRSon'), Qon = clean(Wlead.QRSon(:)); end

    if isfield(Wlead,'QRS_off'), Qoff = clean(Wlead.QRS_off(:));
    elseif isfield(Wlead,'QRSoff'), Qoff = clean(Wlead.QRSoff(:)); end

    if isfield(Wlead,'T_on'), Ton = clean(Wlead.T_on(:));
    elseif isfield(Wlead,'Ton'), Ton = clean(Wlead.Ton(:)); end

    if isfield(Wlead,'T'), Tpk = clean(Wlead.T(:));
    elseif isfield(Wlead,'Tpos'), Tpk = clean(Wlead.Tpos(:)); end

    if isfield(Wlead,'T_off'), Toff = clean(Wlead.T_off(:));
    elseif isfield(Wlead,'Toff'), Toff = clean(Wlead.Toff(:)); end
end

function idx = pick_before(arr, ref, maxBack, minBack)
    % pick latest marker before ref in a valid window
    idx = NaN;
    if isempty(arr), return; end
    cand = arr(arr <= ref - minBack & arr >= ref - maxBack);
    if isempty(cand), return; end
    idx = cand(end);
end

function idx = pick_after(arr, ref, maxFwd, minFwd)
    % pick earliest marker after ref in a valid window
    idx = NaN;
    if isempty(arr), return; end
    cand = arr(arr >= ref + minFwd & arr <= ref + maxFwd);
    if isempty(cand), return; end
    idx = cand(1);
end

function biomFile = find_first_marker_biom_file(rootPa, patName, dateName, sep)
    % find biomarker mat file in Markers folder
    biomFile = '';
    markDir = [rootPa sep patName sep dateName sep 'Markers'];
    if ~exist(markDir,'dir'), return; end

    d = dir([markDir sep '*_Beat3_Biomarkers.mat']);
    if isempty(d), return; end

    biomFile = [markDir sep d(1).name];
end

function dateNameOut = find_first_usable_folder_in_year(rootPa, patName, dateNames, dateDT, targetYear, coreFields, minValidLeads, sep)
    % inside a target year, find the first usable date folder
    dateNameOut = '';

    idxYear = find(year(dateDT) == targetYear);
    if isempty(idxYear), return; end

    for ii = 1:numel(idxYear)
        dn = dateNames{idxYear(ii)};
        biomFile = find_first_marker_biom_file(rootPa, patName, dn, sep);
        if isempty(biomFile), continue; end

        S = load(biomFile);
        if ~isfield(S,'biomarkers'), continue; end
        B = S.biomarkers;

        okFields = true;
        for cf = 1:numel(coreFields)
            if ~isfield(B, coreFields{cf}), okFields = false; break; end
        end
        if ~okFields, continue; end

        valid = true(1,12);
        for cf = 1:numel(coreFields)
            v = B.(coreFields{cf});
            valid = valid & ~isnan(v);
        end

        if sum(valid) >= minValidLeads
            dateNameOut = dn;
            return;
        end
    end
end

function [h,p,note] = safe_lillie(x, alphaLevel)
    % wrapper around lillietest
    h = NaN; p = NaN; note = 'OK';

    x = x(:);
    x = x(~isnan(x));

    if numel(x) < 4
        note = 'TooFewSamples';
        return;
    end

    try
        [h,p] = lillietest(x,'Alpha',alphaLevel);
    catch ME
        note = ['lillietestFailed: ' ME.message];
        h = NaN; p = NaN;
    end
end

function [Y1_med, Y2_med, D_med] = init_stat_structs(fieldList)
    % initialize empty structs for patient-level medians and differences
    Y1_med = struct();
    Y2_med = struct();
    D_med  = struct();

    for f = 1:numel(fieldList)
        Y1_med.(fieldList{f}) = [];
        Y2_med.(fieldList{f}) = [];
        D_med.(fieldList{f})  = [];
    end
end

function [TestT, forceWilcoxon] = build_global_test_table(fieldList, Y1_med, Y2_med, D_med, alphaLevel)
    % build final stats table for one comparison
    %
    % PASS 1:
    % check normality of all biomarker differences
    %
    % PASS 2:
    % if any biomarker is non-normal / fails -> Wilcoxon for all
    % otherwise -> paired t-test for all

    lillieInfo = struct();
    forceWilcoxon = false;

    % pass 1: normality / failure check
    for f = 1:numel(fieldList)
        field = fieldList{f};

        d  = D_med.(field);
        y1 = Y1_med.(field);
        y2 = Y2_med.(field);

        if isempty(d) || isempty(y1) || isempty(y2), continue; end
        if ~(numel(d)==numel(y1) && numel(d)==numel(y2)), continue; end

        [hL,pL,noteL] = safe_lillie(d, alphaLevel);

        lillieInfo.(field).h = hL;
        lillieInfo.(field).p = pL;
        lillieInfo.(field).note = noteL;

        if contains(noteL, 'lillietestFailed', 'IgnoreCase', true) || (isfinite(hL) && hL == 1)
            forceWilcoxon = true;
        end
    end

    % pass 2: paired tests
    testRows = cell(0,12);

    for f = 1:numel(fieldList)
        field = fieldList{f};

        d  = D_med.(field);
        y1 = Y1_med.(field);
        y2 = Y2_med.(field);

        if isempty(d) || isempty(y1) || isempty(y2), continue; end
        if ~(numel(d)==numel(y1) && numel(d)==numel(y2)), continue; end

        if isfield(lillieInfo, field)
            hL = lillieInfo.(field).h;
            pL = lillieInfo.(field).p;
            noteL = lillieInfo.(field).note;
        else
            [hL,pL,noteL] = safe_lillie(d, alphaLevel);
        end

        meanDiff   = mean(d,'omitnan');
        medianDiff = median(d,'omitnan');
        meanY1     = mean(y1,'omitnan');
        meanY2     = mean(y2,'omitnan');

        if forceWilcoxon
            testUsed = 'signrank (global)';
            try
                pT = signrank(y2, y1);
                hT = double(pT < alphaLevel);
                testNote = ['GlobalWilcoxon: ' noteL];
            catch ME
                pT = NaN; hT = NaN;
                testNote = ['signrankFailed: ' ME.message];
            end
        else
            testUsed = 'ttest';
            try
                [hT,pT] = ttest(y2, y1, 'Alpha', alphaLevel);
                testNote = noteL;
            catch ME
                pT = NaN; hT = NaN;
                testNote = ['ttestFailed: ' ME.message];
            end
        end

        testRows(end+1,:) = { ...
            field, numel(d), hL, pL, testNote, ...
            testUsed, hT, pT, meanDiff, medianDiff, meanY1, meanY2}; %#ok<AGROW>
    end

    testVarNames = {'Biomarker','N_patients','Lillie_h_Diff','Lillie_p_Diff','NormalityNote', ...
                    'PairedTest','Test_h','Test_p','MeanDiff_Y2minusY1','MedianDiff_Y2minusY1','MeanY1','MeanY2'};

    if isempty(testRows)
        TestT = cell2table(cell(0,numel(testVarNames)), 'VariableNames', testVarNames);
    else
        TestT = cell2table(testRows, 'VariableNames', testVarNames);
    end
end

function ylabel_for_field(field)
    % set y-axis label depending on biomarker type
    if contains(field,'_ms')
        ylabel('\Delta [ms]');
    elseif contains(field,'AreaQRS_abs')
        ylabel('\Delta [uV*s approx]');
    elseif contains(field,'_mV') || contains(field,'ampl') || contains(field,'Peak')
        ylabel('\Delta [mV]');
    else
        ylabel(['\Delta ' field],'Interpreter','none');
    end
end

function out = ternary_str(cond, a, b)
    % if/else helper
    if cond
        out = a;
    else
        out = b;
    end
end