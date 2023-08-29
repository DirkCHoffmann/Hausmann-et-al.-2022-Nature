function [main_analysis, trigger, hubs] = hubcells (C, S, cfile, im, threshcorr, threshtriggercells, output, main_analysis, scaling, tracesleft, threshhubcells)
posShifts = zeros(1,size(S,1));
correlations = zeros(1,size(S,1));

for i = 1 : size(S,1)
    shifts_of_roi = S(:,i);
    for ii = 1 : size(S,1)
        if C(ii,i) >= threshcorr
            correlations(i) = correlations(i)+1;
            if shifts_of_roi(ii)>0
                posShifts(i) = posShifts(i)+1;
            end
        end
    end
end

n=0;
for i=1 : size(posShifts,2)
    if posShifts(i) >= threshtriggercells
        n=n+1;
        trigger(i) = 1;
        xt(n)= cfile(1,(i*2)-1)/scaling;
        yt(n)= cfile(1,(i*2))/scaling;
    else
        trigger(i) = 0;
    end
end

cellnumber = size(cfile,2)/2;
x1=zeros(cellnumber,1);
y1=zeros(cellnumber,1);
for i=1 : cellnumber
    x1(i)= cfile(1,(i*2)-1)/scaling;
    y1(i)= cfile(1,(i*2))/scaling;
end

figure
imshow(im); %öffnet das Bild
hold on %sorgt dafür dass es bearbeitet wird und kein neuer Plot erstellt wird
b=plot(x1,y1,'g.');
set(b,'MarkerSize',5);
if exist('xt', 'var')
    b=plot(xt,yt,'b.');
    set(b,'MarkerSize',15);
end
drawnow
savefig(output + "_triggercells.fig")

n=0;
for i=1 : size(correlations,2)
    if correlations(i) >= threshhubcells
        n=n+1;
        hubs(i) = 1;
        xh(n)= cfile(1,(i*2)-1)/scaling;
        yh(n)= cfile(1,(i*2))/scaling;
    else
        hubs(i) = 0;
    end
end

figure
imshow(im); %öffnet das Bild
hold on %sorgt dafür dass es bearbeitet wird und kein neuer Plot erstellt wird
b=plot(x1,y1,'g.');
set(b,'MarkerSize',5);
if exist('xh', 'var')
    b=plot(xh,yh,'b.');
    set(b,'MarkerSize',15);
end
drawnow
savefig(output + "_hubcells.fig")


main_analysis(14,1) = "trigger cells";
main_analysis(14,2) = sum(trigger(trigger==1));
main_analysis(15,1) = "proportion of trigger cells (#14/#6)";
if tracesleft == 0
    main_analysis(15,2) = "NaN";
else
    main_analysis(15,2) = sum(trigger(trigger==1)) / tracesleft;
end

main_analysis(16,1) = "hub cells";
main_analysis(16,2) = sum(hubs(hubs==1));
main_analysis(17,1) = "proportion of hub cells (#16/#6)";
if tracesleft == 0
    main_analysis(17,2) = "NaN";
else
    main_analysis(17,2) = sum(hubs(hubs==1)) / tracesleft;
end

trigger = sparse(trigger);
writematrix (trigger, output + "_triggercells.txt");
hubs = sparse(hubs);
writematrix (hubs, output + "_hubcells.txt");

end