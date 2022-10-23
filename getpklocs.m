function [pklocs_long, main_analysis, amplitudes] = getpklocs(ifile, minamplitude, maxwidth, seconds_per_frame, output, show_traces, firsttrace, lasttrace, im, cfile, scaling)

frames = size(ifile,1);
rois = size(ifile,2);
pklocs_long = zeros(size(ifile,1),size(ifile,2));

framespermin = 60/seconds_per_frame;


%Zieht von jedem Rois durch den Durchschnitt aller Rois ab

framemean=zeros(1,size(ifile,2));
for frame=1 : size(ifile,1)
    framemean(1,frame) = mean(ifile(frame,:)); %Durchschnitt der Zeile = Durchschnitt des Frames
    for roi=1 : size(ifile,2)
        ifile(frame,roi) = ifile(frame,roi) - framemean(1,frame); %Intensität der Zelle - Durchschnitt des Frames
    end
end

ifile = gaussfilt(ifile, seconds_per_frame);

%ifile_vector = ifile(:);
%deviation = std(ifile_vector);

for trace = 1:size(ifile,2)
    
    %jetzt wird 'findpeaks' ausgeführt
    %wir bekommen die Höhe, die Position, die Breite und die Amplitude der
    %Peaks (height,locs, width, amp):
    
    [height,locs, width, amp]   = findpeaks(ifile(:,trace),'MinPeakProminence',(minamplitude), 'MaxPeakWidth', (maxwidth/seconds_per_frame));
    PeakInfo(trace).Pklocs      = locs;
    PeakInfo(trace).Pkheight    = ifile(locs,trace);
    
    amplitudes(locs,trace)      = amp;
    widths(locs,trace)          = width;
    
end

pklocs = amplitudes > 0;
s=size(pklocs,1);
if s<frames
    pklocs_long(1:s,:) = pklocs;
end

writematrix (pklocs, output + "_pklocs.txt");

%Anzahl der peaks:
peaksum = sum(sum(amplitudes > 0));
peaks_per_10min = peaksum / frames * framespermin*10;
peaks_norm = peaks_per_10min / rois * 100;

%durchschnittliche Amplitude aller Peaks aus allen Rois:
peakamplitudes = reshape(amplitudes(pklocs),1,[]).';
ampmean = mean(peakamplitudes);
writematrix (peakamplitudes, output + "_peakamplitudes.txt");


%durchschnittliche Peak Breite (=Dauer):
peakwidths = reshape(widths(pklocs),1,[]);
meanwidth = mean(peakwidths, 'all');
peakwidths = peakwidths.' * seconds_per_frame;
writematrix (peakwidths, output + "_peakwidths.txt");

%sprintf('Frames = %.0f \n Rois = %.0f',frames, rois)
%sprintf('Main analysis: \n peaks / 10min & 100 cells = %.1f \n Mean amplitude: %.1f \n Mean width of amplitudes = %.1f seconds', peaks_norm, ampmean, meanwidth)

main_analysis(1,1) = "Frames";
main_analysis(1,2) = frames;
main_analysis(2,1) = "Rois";
main_analysis(2,2) = rois;
main_analysis(3,1) = "peaks / 10min & 100 cells";
main_analysis(3,2) = peaks_norm;
main_analysis(4,1) = "Mean peak amplitude";
main_analysis(4,2) = ampmean;
main_analysis(5,1) = "Mean peak width (seconds)";
main_analysis(5,2) = meanwidth;

if show_traces > 0
    for trace = firsttrace : lasttrace
        figure()
        plot(ifile(:,trace))
        ylim([-100 100]);
        hold on
        plot(PeakInfo(trace).Pklocs,PeakInfo(trace).Pkheight, 'o')
        hold off
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', round(xt*seconds_per_frame))
        drawnow
    end
    
    cellnumber = size(cfile,2)/2;
    x=zeros(cellnumber,1);
    y=zeros(cellnumber,1);
    for i=firsttrace : lasttrace
        x(i)= cfile(1,(i*2)-1)/scaling;
        y(i)= cfile(1,(i*2))/scaling;
        l(i)= i;
    end
    
    figure
    imshow(im); %öffnet das Bild
    hold on %sorgt dafür dass es bearbeitet wird und kein neuer Plot erstellt wird
    b=plot(x,y,'b.');
    set(b,'MarkerSize',15);
    tx = 10; ty = -10;  % displacement so the text does not overlay the data points
    for i=1 : size(l,2)
        t = text(x(i)+tx, y(i)+ty, num2str(round(l(i))));
        set(t, 'Color',[1, 0 ,0], 'FontSize', 15)
    end
    drawnow
    
end


end

function ifile=gaussfilt(ifile, seconds_per_frame)
for g=1 : size(ifile,2)
    ifile(:,g) = imgaussfilt(ifile(:,g),10/seconds_per_frame); %% default: "10"
end
end