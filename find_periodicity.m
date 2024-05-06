function [main_analysis,periodicity] = find_periodicity (pklocs, cfile, im, thresh, output, main_analysis, scaling, do_triggercells, do_correlations, trigger, hubs, seconds_per_frame, S, C, threshcorr, pklim, Dir, amplitudes)
%Kürzen von pklocs auf XX Minütige Sequenzen
lengthseq = 60;    %Länge der Sequenzen die korreliert werden in Minuten default: 60
steps = 5;  %die Zeit, um die die Sequenzen verschoben werden default: 5
numofseq = ceil(((size(pklocs,1)*seconds_per_frame/60)-lengthseq) / steps); %gibt die Anzahl der Seqenzen an, die angeschaut werden sollen: (Länge des Videos - die Länge der letzen Sequenz)/die Zeit, um die die Sequenzen verschoben werden
if numofseq < 1
    numofseq = 1;
end

for a=1 : numofseq
    start = floor(((a-1)*(steps*60/seconds_per_frame)))+1;
    lim = ceil(start-1+(lengthseq*60/seconds_per_frame));
    if a == numofseq
        lim = size(pklocs,1);
    end
    di=lim-start;
    seq_pklocs(a,1:di+1,:) = pklocs(start:lim,:);
end

meanperiod = NaN(1,size(pklocs,2));
deviation = NaN(1,size(pklocs,2));
meanperiodsteps = NaN(numofseq,size(pklocs,2));
deviationsteps = NaN(numofseq,size(pklocs,2));
cellsabove3peaks = 0;

for roi = 1 : size(pklocs,2)
    
    if sum(pklocs(:,roi)) >= pklim
        
        cellsabove3peaks = cellsabove3peaks + 1;
        
        for step=1 : numofseq %Diese Schleife (step) wiederholt die Korrelationen für jede Sequenz (seq) und wählt am ende die niedrigste std aus aus
            
            pkpos = find(seq_pklocs(step,:,roi));
            if sum(seq_pklocs(step,:,roi)) > 3
                
                period = zeros(1,size(pkpos,1)-1);
                for i = 1 : size(pkpos,2)-1
                    period(i) = pkpos(i+1) - pkpos(i);
                end
                
                meanperiodsteps(step,roi) = mean(period);
                deviationsteps(step,roi) = std(period);
                
            else
                meanperiodsteps(step,roi) = NaN;
                deviationsteps(step,roi) = NaN;
            end
            
            clear pkpos period
            
        end
        
        [deviation(roi),loc] = min(deviationsteps(:,roi));
        meanperiod(roi) = meanperiodsteps(loc,roi);
        
        
    else
        meanperiod(roi) = NaN;
        deviation(roi) = NaN;
        
    end
end

if do_correlations > 0
    % hier werden die Anzahl der Korrelationen der periodischen und
    % nicht-periodischn Zellen gezählt
    corrall = NaN(1,size(S,1));
    posShiftsall = NaN(1,size(S,1));
    posShiftsperio = NaN(1,size(S,1)); %Anzahl der Zellen die getriggert werden (aller periodischen Zellen)
    corrperio = NaN(1,size(S,1)); %Anzahl der verbundenen Zellen (aller periodischen Zellen)
    posShiftsNOTperio = NaN(1,size(S,1)); %Anzahl der Zellen die getriggert werden (aller NICHT periodischen Zellen)
    corrNOTperio = NaN(1,size(S,1)); %Anzahl der verbundenen Zellen (aller NICHT periodischen Zellen)
    dir_all = NaN(size(S,1)); %Direktionalitäten aller periodischen Zellen
    dir_perio = NaN(size(S,1)); %Direktionalitäten aller periodischen Zellen
    dir_NOTperio = NaN(size(S,1)); %Direktionalitäten aller NICHT periodischen Zellen
    
    for i = 1 : size(S,1)
        if sum(pklocs(:,i)) > 3 % Die Zelle muss mindestens 4x Blinken, um als periodische oder nicht-periodische Zelle Eingestuft zu werde
            
            shifts_of_roi = S(i,:);
            if deviation(i) <= thresh/seconds_per_frame % wenn es eine periodische Zelle ist:
                if isnan(corrall(i))
                    corrall(i)=0;
                    posShiftsall(i) = 0;
                    corrperio(i) = 0;
                    posShiftsperio(i) = 0;
                end
                
                for ii = 1 : size(S,1)
                    if C(ii,i) >= threshcorr
                        corrperio(i) = corrperio(i)+1;
                        corrall(i) = corrall(i)+1;
                        if shifts_of_roi(ii)>0
                            posShiftsperio(i) = posShiftsperio(i)+1;
                            posShiftsall(i) = posShiftsall(i)+1;
                        end
                    end
                end
                dir_perio(i,:) = Dir(i,:);
                dir_all(i,:) = Dir(i,:);
                
            else % wenn es eine NICHT periodische Zelle ist:
                if isnan(corrall(i))
                    corrall(i)=0;
                    posShiftsall(i) = 0;
                    corrNOTperio(i) = 0;
                    posShiftsNOTperio(i) = 0;
                end
                for ii = 1 : size(S,1)
                    if C(ii,i) >= threshcorr
                        corrNOTperio(i) = corrNOTperio(i)+1;
                        corrall(i) = corrall(i)+1;
                        if shifts_of_roi(ii)>0
                            posShiftsNOTperio(i) = posShiftsNOTperio(i)+1;
                            posShiftsall(i) = posShiftsall(i)+1;
                        end
                    end
                end
                dir_NOTperio(i,:) = Dir(i,:);
                dir_all(i,:) = Dir(i,:);
                
            end
        end
    end
    
    %dir_perio = (posShiftsperio + posShiftsperio - corrperio)./corrperio;
    %dir_NOTperio = (posShiftsNOTperio + posShiftsNOTperio - corrNOTperio)./corrNOTperio;
    
    
    main_analysis(24,1) = "Mean number of correlating cells per periodic cell";
    main_analysis(24,2) = mean(corrperio,'omitnan');
    main_analysis(25,1) = "Mean number of correlating cells per NOT-periodic cell that have at least 4 peaks";
    main_analysis(25,2) = mean(corrNOTperio,'omitnan');
    main_analysis(26,1) = "Mean number of triggered cells per periodic cell";
    main_analysis(26,2) = mean(posShiftsperio,'omitnan');
    main_analysis(27,1) = "Mean number of triggered cells per NOT-periodic cell that have at least 4 peaks";
    main_analysis(27,2) = mean(posShiftsNOTperio,'omitnan');

    dir_all = dir_all(:);
    dir_all = dir_all(~isnan(dir_all));
    dir_perio = dir_perio(:);
    dir_perio = dir_perio(~isnan(dir_perio));
    dir_NOTperio = dir_NOTperio(:);
    dir_NOTperio = dir_NOTperio(~isnan(dir_NOTperio));
    
    corrall = corrall(~isnan(corrall));
    posShiftsall = posShiftsall(~isnan(posShiftsall));
    corrperio = corrperio(~isnan(corrperio));
    corrNOTperio = corrNOTperio(~isnan(corrNOTperio));
    posShiftsperio = posShiftsperio(~isnan(posShiftsperio));
    posShiftsNOTperio = posShiftsNOTperio(~isnan(posShiftsNOTperio));
    net_all = (2*posShiftsall)-corrall;
    net_perio = (2*posShiftsperio)-corrperio;
    net_NOTperio = (2*posShiftsNOTperio)-corrNOTperio;
    
    if size(dir_all,1) > size(net_all,2)
        mainmatrix = NaN(size(dir_all,1),10);
    else
        mainmatrix = NaN(size(net_all,2),10);
    end
    mainmatrix(1:size(corrall,2),1) = corrall.';
    mainmatrix(1:size(corrperio,2),2) = corrperio.';
    mainmatrix(1:size(corrNOTperio,2),3) = corrNOTperio.';
    mainmatrix(1:size(net_all,2),4) = net_all.';
    mainmatrix(1:size(net_perio,2),5) = net_perio.';
    mainmatrix(1:size(net_NOTperio,2),6) = net_NOTperio.';
    mainmatrix(1:size(dir_all),8) = dir_all;
    mainmatrix(1:size(dir_perio),9) = dir_perio;
    mainmatrix(1:size(dir_NOTperio),10) = dir_NOTperio;
    
    writematrix (mainmatrix, output + "_corrall-corrperio-corrNOTperio-netall-net_perio-net_NOTperio-dirall-dirperio-dirnotperio.txt");
end

cellnumber = size(cfile,2)/2;
periodicity = zeros(1,cellnumber);
n=0;
m=0;
amplitudes_perio = NaN;
amplitudes_non_perio = NaN;

for i=1 : cellnumber
    if deviation(i) <= thresh/seconds_per_frame
        n=n+1;
        periodicity(i) = 1;
        periodsofperiodiccells(n) = meanperiod(i)*seconds_per_frame;
        x(n) = cfile(1,(i*2)-1)/scaling;
        y(n) = cfile(1,(i*2))/scaling;
        p(n) = meanperiod(i)*seconds_per_frame;
        
        activityperiodic(n) = sum(pklocs(:,i));
        amplitudes_perio(n) = mean(nonzeros(amplitudes(:,i)));
        
    elseif deviation(i) > thresh/seconds_per_frame
        m=m+1;
        activitynotperiodic(m) = sum(pklocs(:,i));
        
        amplitudes_non_perio(m) = mean(nonzeros(amplitudes(:,i)));

    end
end

writematrix (amplitudes_perio.', output + "_peakamplitudes_perio.txt");
writematrix (amplitudes_non_perio.', output + "_peakamplitudes_non-perio.txt");

if exist('periodicity','var')
else
    periodicity = 0;
end

x1=zeros(cellnumber,1);
y1=zeros(cellnumber,1);
for i=1 : cellnumber
    x1(i)= cfile(1,(i*2)-1)/scaling;
    y1(i)= cfile(1,(i*2))/scaling;
end

figure('Name','Periodic cells');
imshow(im); %öffnet das Bild
hold on %sorgt dafür dass es bearbeitet wird und kein neuer Plot erstellt wird
b=plot(x1,y1,'g.');
set(b,'MarkerSize',5);
if exist('x', 'var')
    b=plot(x,y,'b.');
    set(b,'MarkerSize',15);
    tx = 10; ty = -10;  % displacement so the text does not overlay the data points
    for i=1 : size(p,2)
        t = text(x(i)+tx, y(i)+ty, num2str(round(p(i))));
        set(t, 'Color',[1, 0 ,0], 'FontSize', 15)
    end
end
hold off
drawnow
savefig(output + "_periodicity.fig")


main_analysis(9,1) = "cells with at least 4 peaks";
main_analysis(9,2) = cellsabove3peaks;
main_analysis(10,1) = "periodic cells";
main_analysis(10,2) = sum(periodicity(periodicity==1));
main_analysis(11,1) = "mean of all periodic periods";
main_analysis(12,1) = "std of all periodic periods";
if exist('periodsofperiodiccells','var')
    main_analysis(11,2) = mean(periodsofperiodiccells);
    main_analysis(12,2) = std(periodsofperiodiccells);
else
    main_analysis(11,2) = "NaN";
    main_analysis(12,2) = "NaN";
    periodsofperiodiccells = NaN;
end
writematrix (periodsofperiodiccells.', output + "_periodsofperiodiccells.txt");


main_analysis(13,1) = "proportion of periodic cells (#10/#9)";
if cellsabove3peaks == 0
    main_analysis(13,2) = "NaN";
else
    main_analysis(13,2) = sum(periodicity(periodicity==1)) / cellsabove3peaks;
end

sparse_periodicity = sparse(periodicity);
sparse_meanperiod = sparse(meanperiod(~isnan(meanperiod)));
sparse_deviation = sparse(deviation(~isnan(deviation)));
combined = sparse_meanperiod*seconds_per_frame;
combined(2,:) = sparse_deviation*seconds_per_frame;
TF = isempty(combined);
if TF == 1
    combined = 0;
end
writematrix (sparse_periodicity, output + "_periodicity.txt");
writematrix (combined, output + "_mean_std_periods.txt");


if do_triggercells > 0 && do_correlations > 0
    trigger = full(trigger);
    hubs = full(hubs);
    
    overlaptri = sum(periodicity(trigger==1)); %zählt Zellen die sowohl trigger als auch periodisch sind
    overlaphub = sum(periodicity(hubs==1)); %zählt Zellen die sowohl hubs als auch periodisch sind
    
    main_analysis(18,1) = "periodic cells in trigger cells";
    main_analysis(18,2) = overlaptri / sum(trigger(trigger==1));
    main_analysis(19,1) = "trigger cells in periodic cells";
    main_analysis(19,2) = overlaptri / sum(periodicity(periodicity==1));

    
    main_analysis(20,1) = "periodic cells in hub cells";
    main_analysis(20,2) = overlaphub / sum(hubs(hubs==1));
    main_analysis(21,1) = "hub cells in  periodic cells";
    main_analysis(21,2) = overlaphub / sum(periodicity(periodicity==1));

end

if exist('activityperiodic','var')
    meanactperio = mean(activityperiodic(activityperiodic>3)); %die durchschnittliche Anzahl der Peaks aller periodischen Zellen
else
    meanactperio = NaN;
    activityperiodic = NaN;
end

if exist('activitynotperiodic','var')
    meanactnotperio = mean(activitynotperiodic(activitynotperiodic>3)); %die durchschnittliche Anzahl der Peaks aller NICHT-periodischen Zellen mit mind. 4 Peaks
else
    meanactnotperio = NaN;
    activitynotperiodic = NaN;
end

main_analysis(22,1) = "Mean number of peaks of periodic cells";
main_analysis(22,2) = meanactperio;
main_analysis(23,1) = "Mean number of peaks of NOT-periodic cells that have at least 4 peaks";
main_analysis(23,2) = meanactnotperio;

perioactivity = NaN(size(pklocs,2),2);
perioactivity(1:size(activityperiodic,2),1) = activityperiodic.';
perioactivity(1:size(activitynotperiodic,2),2) = activitynotperiodic.';
writematrix (perioactivity, output + "_perioactivity-perio-notperio.txt");


% Trigger coefficiant and Cluster coefficiant
if do_triggercells > 0 && do_correlations > 0
    
    triggercoeff = NaN(size(pklocs,2),1);
    clustercoeff = NaN(size(pklocs,2),1);
    
    allcorrelations = length(find(C>=threshcorr))/2;
    alltriggers = allcorrelations/2;
    
    mean_triggering = alltriggers / cellsabove3peaks;
    mean_clustering = allcorrelations / cellsabove3peaks;
    
    for i=1 : size(pklocs,2)
        if sum(pklocs(:,i))>=pklim
            triggercoeff(i,1) = length(find(C(i,S(i,:)>0)>=threshcorr)) / (2*mean_triggering);
            clustercoeff(i,1) = length(find(C(i,:)>=threshcorr)) / (2*mean_clustering);
        end
    end
    
    triggercoeff_perio = triggercoeff(find(deviation <= thresh/seconds_per_frame),1);
    triggercoeff_notperio = triggercoeff(find(deviation > thresh/seconds_per_frame),1);
    clustercoeff_perio = clustercoeff(find(deviation <= thresh/seconds_per_frame),1);
    clustercoeff_notperio = clustercoeff(find(deviation > thresh/seconds_per_frame),1);
    
    perio_coeff = NaN(size(pklocs,2),4);
    perio_coeff(1:size(clustercoeff_perio,1),1) = clustercoeff_perio;
    perio_coeff(1:size(clustercoeff_notperio,1),2) = clustercoeff_notperio;
    perio_coeff(1:size(triggercoeff_perio,1),3) = triggercoeff_perio;
    perio_coeff(1:size(triggercoeff_notperio,1),4) = triggercoeff_notperio;
    
    
    writematrix (perio_coeff, output + "_perio_coeff.txt");
end
end