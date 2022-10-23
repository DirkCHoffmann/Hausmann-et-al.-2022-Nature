function lines_on_image(im ,cfile, C, threshcorr, output, scaling, numbering)

cellnumber = size(cfile,2)/2;
x=zeros(cellnumber,1);
y=zeros(cellnumber,1);
for i=1 : cellnumber
    x(i)= cfile(1,(i*2)-1)/scaling;
    y(i)= cfile(1,(i*2))/scaling;
end

ColorMap = jet(100); %Verteilt die Farbscala "jet" auf die Zahlen 0 bis 100

figure
imshow(im); %öffnet das Bild
hold on %sorgt dafür dass es bearbeitet wird und kein neuer Plot erstellt wird
a=plot(x,y,'r.');
set(a,'MarkerSize',10);
%drawnow
for i=1 : cellnumber
    for ii=1 : cellnumber
        if C(ii,i) >= threshcorr
            %i
            %ii
            LineColor = ColorMap(fix(C(ii,i)*100), :); %Die Linecolor entspricht der Korrelation*100 auf der jet scala
            b=line([x(i) x(ii)], [y(i) y(ii)], 'Color', LineColor);
            set(b,'LineWidth',2);
            
        end
    end
    if numbering > 0
    tx = 5; ty = -5;  % displacement so the text does not overlay the data points
    t = text(x(i)+tx, y(i)+ty, num2str(i));
    set(t, 'Color',[1, 0 ,0], 'FontSize', 6)
    end
end

drawnow
output_figure = sprintf('%s_lines on image.fig', output);
savefig(output_figure);

end

