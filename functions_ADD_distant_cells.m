
function [D, C, S, corrleft, tracesleft, main_analysis, Dir] = functions_ADD_distant_cells(cfile, ifile, pklocs, minspeed, maxspeed, dcutoff, pklim, output, seconds_per_frame, main_analysis, threshcorr, do_null)
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% festlegen welche Zellen nicht berücksichtig werden sollen
not=zeros(size(ifile,2),1);
for z=1 : size(ifile,2)
    if sum(pklocs(:,z)) < pklim % die Zellen mit weniger peaks als peaklim sollen nicht berücksichtigt werden
        not(z)=1;
    end
end

alltraces = size(ifile,2);
below_pklim = sum(not) / alltraces * 100;
tracesleft = alltraces-sum(not);
allcorr = ((tracesleft * tracesleft) - tracesleft ) / 2;

main_analysis(6,1) = "Rois with number of peaks not below pklim";
main_analysis(6,2) = tracesleft;

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% DISTANZ: %berechnet die Distanzen zwischen jedem Punkt als Matrix
cellnumber = size(cfile,2)/2;
D=NaN(size(cellnumber));
%cfile wird in X(i) und Y(ii) werte übertragen
corrleft = 0;
for i=1 : cellnumber
    for ii=1 : cellnumber
        xi= cfile(1,(i*2)-1);
        yi= cfile(1,(i*2));
        xii= cfile(1,(ii*2)-1);
        yii= cfile(1,(ii*2));
        
        D(i,ii) = sqrt(((xi-xii)^2)+((yi-yii)^2)); %aus den X- und Y-Werten der beiden Zellen wird die Distanz berechnet
        
        if not(i) == 0 && not(ii)==0 && D(i,ii) > dcutoff &&  i~=ii %die anzahl aller Zellpaare unter dcutoff, bei der beide Zellen mindestens 4 peaks haben
            corrleft = corrleft + 1;
        end
    end
end
corrbelowdcutoff1 = size(D((0<D)&(D>dcutoff))); %die anzahl aller Zellpaare unter dcutoff
corrbelowdcutoff2 = corrbelowdcutoff1(1)/2;
corrleft = corrleft/2;

main_analysis(7,1) = "all possible correlations below dcutoff";
main_analysis(7,2) = corrbelowdcutoff2;
main_analysis(8,1) = "corrleft (all possible correlations below dcutoff between cells with number of peaks not below pklim)";
main_analysis(8,2) = corrleft;

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% berechnet die maximale Korrelation durch Verschiebung

%Kürzen von ifile auf 10 Minütige Sequenzen
length = 10;    %Länge der Sequenzen die korreliert werden in Minuten
steps = 5;  %die Zeit, um die die Sequenzen verschoben werden
numofseq = ceil(((size(ifile,1)*seconds_per_frame/60)-length) / steps); %gibt die Anzahl der Seqenzen an, die korreliert werden sollen: (Länge des Videos - die Länge der letzen Sequenz)/die Zeit, um die die Sequenzen verschoben werden
if numofseq<1
    numofseq=1;
end
fprintf('\n\nNumber of Sequences: %.0f \nLength of Sequence: %.0f\n\n', numofseq, length)

for a=1 : numofseq
    start = floor(((a-1)*(steps*60/seconds_per_frame)))+1;
    lim = ceil(start-1+(length*60/seconds_per_frame));
    if a == numofseq
        lim = size(ifile,1);
    end
    di=lim-start;
    seq(a,1:di+1,:) = ifile(start:lim,:);
    seq_pklocs(a,1:di+1,:) = pklocs(start:lim,:);
end

%Erstellen aller benötigten Matizen, damit sie nicht in den folgenen Loops
%immer wieder ihre Größe verändern müssen:
X=zeros(size(ifile,2),size(ifile,2));
Cs=NaN(numofseq,size(ifile,2),size(ifile,2));
Ss=NaN(numofseq,size(ifile,2),size(ifile,2));
Ps=NaN(numofseq,size(ifile,2),size(ifile,2));
C=NaN(size(ifile,2));
S=NaN(size(ifile,2));
P=NaN(size(ifile,2));
Dir=NaN(size(ifile,2));
corr_pos_max=NaN(size(ifile,2));
corr_neg_max=NaN(size(ifile,2));

maxspeed = maxspeed*seconds_per_frame; %Umrechnung in distance/frame
minspeed = minspeed*seconds_per_frame;
maxshift = ceil(max(max(D))/minspeed); % Um die Größe von A und pval schon vor der loop festzuhalten
A=NaN(size(ifile,2),size(ifile,2),(maxshift * 2)+1);
B=NaN(size(ifile,2),size(ifile,2),(maxshift * 2)+1);
pval=zeros(size(ifile,2),size(ifile,2),(maxshift * 2)+1);

twentypercent = round(size(ifile,2)/5);
msg=msgbox(sprintf('Progress: 0 percent'));
%geht mit den ersten beiden for-Schleifen (i & ii) durch die möglichen Zellpaare
for i=1 : size(ifile,2)
    
    %gibt den Fortschritt in 10%-Schritten an
    if rem(i, twentypercent)==0
        delete(msg);
        msg=msgbox(sprintf('Progress: %.0f percent', i/twentypercent*20));
    end
    
    for ii=1 : size(ifile,2)
        
        if i ~= ii  && X(i,ii)==0 && D(ii,i)>dcutoff && not(i) == 0 && not(ii)==0
            %minshift= fix(D(ii,i)/maxspeed);
            %maxshift= ceil(D(ii,i)/minspeed);
            minshift= fix(50/maxspeed);
            maxshift= ceil(50/minspeed);
            
            for step=1 : numofseq %Diese Schleife (step) wiederholt die Korrelationen für jede Sequenz (seq) und wählt am ende die höchste Korrelation aus
                %die dritte for-Schleife (iii) shiftet:
                %Sie beginnt bei keinem shift und shiftet dann abwechselnd nach
                %vorne und nach hinten in steigender Weite (Damit der erste
                %höchste Wert (wenn es mehrere gleich hohe höchste
                %Korrelationen gibt) näher an 0 shift liegt.
                
                step_of_other_vector = step;
                
                if do_null == 3 % vom zweiten trace wird ein zufällig ausgewählter anderer step ausgewählt
                    for kk = 0 : 1000
                        step_of_other_vector = floor(1+rand*numofseq);
                        if step_of_other_vector > numofseq
                            step_of_other_vector = numofseq;
                        end
                        if step_of_other_vector ~= step
                            break
                        end
                    end
                end
                
                nn=0;
                mm=0;
                if sum(seq_pklocs(step,:,i)) >= pklim  && sum(seq_pklocs(step,:,ii)) >= pklim  % die Zellen mit weniger peaks als pklim sollen nicht berücksichtigt werden
                    
                    for iii=(minshift*2) : (maxshift*2)+1
                        if iii>=2
                            if rem(iii, 2) == 0
                                shift=iii/2;
                                shiftpn = 1; % 1 = positiver shift
                            elseif rem(iii, 2) == 1
                                shift=-(iii-1)/2;
                                shiftpn = 0; % 0 = negativer shift
                            end
                            
                            movedvector = circshift(seq(step,:,i),shift);
                            
                            movedvector2 = movedvector.'; %tauscht Reihe mit Spalten, damit corr den ganzen vector am Stück vergleicht
                            othervector = seq(step_of_other_vector,:,ii);
                            othervector2 = othervector.';
                            A(ii,i,iii) = corr(movedvector2, othervector2); %gibt jeweils die Korrelation und ihren P-Wert an
                            B(ii,i,iii) = shift;
                            
                            if shiftpn == 1
                                nn = nn+1;
                                corr_pos(nn,step,ii,i) = A(ii,i,iii);
                                
                            end
                            if shiftpn == 0
                                mm = mm+1;
                                corr_neg(nn,step,ii,i) = A(ii,i,iii);
                            end
                            
                        end
                    end
                    
                    [Cs(step,ii,i), maxiii] = max(A(ii,i,:)); %Die höchste Korrelation wird Cs zugeordnet. Der Ort (iii) dieser Korrelation wird cshift zugeordnet
                    Cs(step,i,ii) = Cs(step,ii,i);
                    Ss(step,i,ii) = B(ii,i,maxiii);
                end
                
            end
            [C(i,ii),loc] = nanmax(Cs(:,ii,i));
            S(i,ii) = Ss(loc,i,ii);
            if C(i,ii) >= threshcorr
                corr_pos_max(i,ii) = nanmax(nanmax(corr_pos(:,:,ii,i)));
                corr_neg_max(i,ii) = nanmax(nanmax(corr_neg(:,:,ii,i)));
                %Dir(i,ii)=(corr_pos_max(i,ii)-corr_neg_max(i,ii));
                Dir(i,ii)=(corr_pos_max(i,ii)-corr_neg_max(i,ii))/C(i,ii)*100;
            end
            
            if isnan(C(i,ii))
                S(i,ii) = NaN;
            end
            
            C(ii,i) = C(i,ii);
            S(ii,i) = - S(i,ii);
            P(ii,i) = P(i,ii);
            Dir(ii,i) = -Dir(i,ii);
            
        end
        
        X(ii,i) = 1;
        
    end
end
delete(msg); %das letzte "Progress-Fenster wird geschlossen
n=0;
XX=zeros(size(C,1));
for i=1 : size(C,1)
    for ii=1 : size(C,2)
        if XX(ii,i)==0 && ii~=i
            n = n+1;
            CvsD(n,1) = C(i,ii);
            CvsD(n,2) = D(i,ii);
            XX(i,ii)=1;
        end
    end
end

writematrix(CvsD, output + "_CvsD.txt");
writematrix(C, output + "_C.txt");

end