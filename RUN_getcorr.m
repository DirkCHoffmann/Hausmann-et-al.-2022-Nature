
%Ordner in dem sich die _mean intensity.txt, die _center.txt und die max
%projection.png Dateien befinden.
cd '/Users/davidhausmann/Desktop/test';

repeat = 0; %wenn repeat = 1 werden keine Dateien geöffnet und keine Korrelationen berechnet sondern die Dateinen und Korrelationen des letzten Durchganges übernommen
if repeat == 0
    clearvars -except repeat;
end
change_corr = 0; % wenn change_corr = 1 wird die veränderte CC.dat Datei eingelesen und die Korrelationen übernommen

% 1 -> Funktion wird durchgeführt.
% 0 -> Funktion wird nicht durchgeführt.

do_null = 0; %DO Not FORGET TO CHANGE THE FOLDER!!
% 0: empirische Daten werden verwendet
% 1: cicrcular shift
% 2: scrambled (alle traces werden jeweils in eine Minute lange Stücke geschnitten und zufällig angeordnet)
% 3: linear shift (best control!)
% 4: weit entfernte Zellen (>400 um) werden als Kontrolle verwendet
dcutoff_distant_cells = 400;

do_celllimit = 0; %Zellzahl einschränken
firstcell = 3;
lastcell = 10;

do_dotplot = 0; %dotplot aller peaks
show_traces = 0; %zeigt die Traces 'firsttrace' bis 'lasttrace' zusammen mit den gefunden peaks
    firsttrace = 30;
    lasttrace = 40;
do_correlations = 1; %findest die Korrelationen zwischen den Zellen
    do_lines_on_image = 1; %Visualisierung der Korrelationen
        numbering = 0; %Nummeriert die Zellen durch
    do_triggercells = 1; %Erkennt die Zellen die am häufigsten zuerst blinken
    do_connectivity = 1; %Gibt alle gefundenen Korrelationen / alle möglichen Korrelationen an
    plot_C_vs_D = 0;
    do_numberofconnections = 1;
    do_networktheory = 1;
    do_networkgraphs = 1;
do_periodicity = 1; %findet die periodisch blinkenen Zellen

seconds_per_frame = 1.56;
scaling = 0.52;

% dotplots, 3Dplots, und Histogramme müssen im Skript unten einzelnd aktiviert werden

%get_corr                                                                                               
spread_s = 20; %Breite der Peaks die korreliert werden in Sekunden                                      
minamplitude = 10; %die Mindestamplitude der Peaks für die Funktion "getpklocs"                         
maxwidth = 115; %die Mindestbreite der Peaks für die Funktion "getpklocs" in Sekunden                  
pklim = 4; %Mindestanzahl der peaks die eine Zelle haben muss um berücksichtigt zu werden.             
minspeed = 5; %in um/s (wird in functions_ADD in um/frames umgerechnet)                                  
maxspeed = 40; %in um/s (wird in functions_ADD in um/frames umgerechnet)                                 
dcutoff = 100; %die maximale Distanz zwischen den Punkten die noch verglichen werden sollen            
threshcorr = 0.47; %Threshold für correlationen die Angezeigt werden sollen                     

%periodicity
thresh_std = 15; %Maximale Standartabweichung der Abstände zwischen den Peaks in Sekunden               default: 15

%triggercells
threshtriggercells = 3; %Wie viele Zellen eine Zelle mindestens erregen muss um als Triggercell zu gelten  
threshhubcells = 5; %mit wieviele Zellen die Zelle mindesens korrelieren muss, um als Hub-Zelle zu gelten   


close_all_windows = 1; %closes all windows after every analyzed time series


% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------

if show_traces == 1
    close_all_windows = 0;
end

if do_null > 0
    do_networkgraphs = 0;
end


% Dateien öffnen und bearbeiten (kürzen)
currentfolder = cd;
listoffiles = dir ('*mean intensities.txt'); %hier kann man das entsprechende Pattern einfügen (* steht für alles)
numberoffiles = length(listoffiles);

for repetition=1 : numberoffiles

    if repeat == 0

        clearvars -except change_corr numbering do_networkgraphs threshhubcells maxwidth spread_s close_all_windows lasttrace firsttrace show_traces scaling seconds_per_frame threshcorr repeat do_connectivity do_triggercells do_numberofconnections plot_C_vs_D threshtriggercells do_connectivity do_lines_on_image do_correlations thresh_std repetition listoffiles numberoffiles framespermin minamplitude pklim minspeed maxspeed dcutoff do_null firstcell lastcell do_celllimit do_dotplot do_periodicity do_networktheory dcutoff_distant_cells;

        path = listoffiles(repetition).name;
        ifile_raw = readmatrix(path);
        ifile = ifile_raw(2:end,2:end);  % bei mean intensities muss noch die erste Spalte und Zeile gelöscht werden!
        %ifile = ifile/256;

        path = path(1:length(path)-21); %löscht "_mean intensities.txt" (löscht die letzten 21 zeichen)
        if do_null == 0
            output = path;
        elseif do_null == 1
            output = path + "_circular shift";
        elseif do_null == 2
            output = path + "_scrambled";
        elseif do_null == 3
            output = path + "_linear shift";
        elseif do_null == 4
            output = path + "_distant cells";
        end
        output_diary = output + "_Report.txt";
        diary(output_diary);
        fprintf('\n\nFile %.0f of %.0f\n\n', repetition, numberoffiles)
        fprintf(path);
        fprintf('\n\nseconds_per_frame: %.3f \nscaling: %.3f ', seconds_per_frame, scaling)

        %einlesen des ..._center.txt file
        center = path + "_center.txt";
        cfile = readmatrix(center);
        cfile = cfile(2,2:end);

        %einlesen des ersten Bildes oder max intensity projection der Zellen
        image = path + "_max projection.png";
        im = imread(image);

        %liefert alle peaks mit einer Mindestamplitude von 'minamplitude' als 20 frames lange Reihe von Einsen
        [pklocs, main_analysis, amplitudes] = getpklocs(ifile, minamplitude, maxwidth, seconds_per_frame, output, show_traces, firsttrace, lasttrace, im, cfile, scaling);

        if do_null == 1
            pklocs = rand_matrix_CircularShift(pklocs); %Für die Kontrolle werden alle traces zufällig cirular geshiftet
        elseif do_null == 2
            pklocs = rand_matrix_Scrambled(pklocs);    %So wurde die Kontrolle bis zur ersten Submission berechnet
        end

        ifile = spread(pklocs,spread_s/seconds_per_frame);

    end

    diary(output_diary);
    %print variables
    fprintf('\n\nminamplitude: %.0f \n\nthreshcorr: %.2f \n\npklim: %.0f \n\nminspeed: %.1f um/s \n\nmaxspeed: %.1f um/s \n\ndistance cutoff: %.0f um', minamplitude, threshcorr, pklim, minspeed, maxspeed, dcutoff)
    fprintf('\n\nthresh_std (Periodicity): %.0f',thresh_std)
    fprintf('\n\nthreshtriggercells: %.0f \n\nthreshhubcells: %.0f', threshtriggercells, threshhubcells)

    if do_celllimit > 0
        %Zellzahl einschränken
        ifile = ifile(:,firstcell:lastcell);
        cfile = cfile(1,(firstcell*2)-1:lastcell*2);
        pklocs = pklocs(:,firstcell:lastcell);
    end

    if do_dotplot > 0
        %dotplot from pkloc
        [frame,cellnumber] = find(ifile);
        b = plot(frame,-cellnumber,'k.');
        set(b,'MarkerSize',5);
        xlabel("Frames");
        ylabel("# of cell");
        xlim([0 max(frame)]);
        ylim([-max(cellnumber) 0]);
        drawnow
        savefig(output + "_pkloc dotplot.fig")
    end

    % Korrelation durch Verschiebung
    if do_correlations > 0
        if repeat == 0
            if do_null == 4
                [D, C, S, corrleft, tracesleft, main_analysis, Dir] = functions_ADD_distant_cells(cfile, ifile, pklocs, minspeed, maxspeed, dcutoff_distant_cells, pklim, output, seconds_per_frame, main_analysis, threshcorr, do_null);

            else
                [D, C, S, corrleft, tracesleft, main_analysis, Dir] = functions_ADD(cfile, ifile, pklocs, minspeed, maxspeed, dcutoff, pklim, output, seconds_per_frame, main_analysis, threshcorr, do_null);
            end
        end

        if change_corr == 1
            sparseCC = path + "_sparseCC.dat";
            sparseCC = load(sparseCC);
            sparseCC = spconvert(sparseCC);
            CC = full(sparseCC);
            CCm=zeros(size(C,1));
            CCm(1:size(CC,1),1:size(CC,2))=CC;
            for i = 1 : size(C,1)
                for ii = 1 : size(C,1)
                    if C(i,ii) >= threshcorr
                        if CCm(i,ii) ~= 1 || CCm(ii,i) ~= 1
                            C(i,ii) = NaN;
                            S(i,ii) = NaN;
                        end
                    end
                end
            end
        end

        if do_lines_on_image > 0
            %Cvector = C(:);
            %threshcorr = prctile(Cvector,90); %threshold für image berechnen
            lines_on_image(im, cfile, C, threshcorr, output, scaling, numbering); %Visualisierung der Korrelationen
        end

        if do_triggercells > 0
            [main_analysis, trigger, hubs] = hubcells (C, S, cfile, im, threshcorr,threshtriggercells, output, main_analysis, scaling, tracesleft, threshhubcells);
        end

        if do_connectivity > 0
            main_analysis = connectivity(C, corrleft, output, main_analysis, threshcorr);
        end

        if plot_C_vs_D > 0
            figure ('Name','Correlation vs. Distance');
            hold on
            a = plot(D,C,'k.');
            set(a,'MarkerSize',10);
            ylabel("correlation");
            xlabel("distance (um)");
            if do_null ~= 4
                xlim([0 dcutoff]);
            end
            ylim([0 1]);
            line([0, dcutoff],[threshcorr threshcorr],'Color','red');
            hold off
            drawnow
            savefig(output + "_corr_vs_distance.fig")
        end

        if do_networktheory > 0
            [main_analysis, CC] = networktheory(C, main_analysis, threshcorr);
            sparseCC = sparse(CC);
            [i,j,val]=find(sparseCC);
            dlmwrite(output + "_sparseCC.dat",[j i val],'delimiter', ' ','newline','pc');
        end

    end

    if do_periodicity > 0
        if do_correlations == 0 || do_triggercells == 0
            trigger = 0;
            hubs = 0;
        end
        if do_correlations == 0
            C=0;
            S=0;
            Dir=0;
        end

        [main_analysis, periodicity] = find_periodicity (pklocs, cfile, im, thresh_std, output, main_analysis, scaling, do_triggercells, do_correlations, trigger, hubs, seconds_per_frame, S, C, threshcorr, pklim, Dir, amplitudes);
    end

    if do_correlations > 0 && do_numberofconnections > 0
        if do_periodicity == 0
            periodicity=0;
        end
        numberofconnections(C,threshcorr, output, periodicity); %erstellt Histogramm der Verbindungen pro Zelle

    end

    if do_networkgraphs > 0 && do_correlations > 0
        %plots the original network and the original network without periodic cells
        %and calculates the larges cluster
        main_analysis = networkgraphs(hubs, periodicity, CC, output, do_periodicity, main_analysis);
    end

    writematrix (main_analysis, output + "_main-analysis.txt");


    diary off

    if close_all_windows > 0
        close all
    end


end
fprintf('\n');
sound(sin(1:600),6000);
%pause(0.2)
%sound(sin(1:600),6000);


function ifilenew=spread(ifile,width)
ifilenew = zeros(size(ifile,1),size(ifile,2));
for g=1 : size(ifile,1)
    for h=1 : size(ifile,2)
        if ifile(g,h) == 1
            ifilenew(g,h)=1;
            for j=1 : width/2
                if g+j <= size(ifile,1)
                    ifilenew(g+j,h) = 1;
                end
                if g-j > 0
                    ifilenew(g-j,h) = 1;
                end
            end
        end
    end
end
end

