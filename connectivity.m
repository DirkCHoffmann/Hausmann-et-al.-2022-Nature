function main_analysis = connectivity(C, corrleft, output, main_analysis, threshcorr)
conn = zeros(11,2);

for i=0 : 10
    cutoff = i/10;
    conn(i+1,1) = cutoff;
    v=size(C(C>=cutoff))/2;
    if v(1)>0
        conn(i+1,2) = v(1) / corrleft;
    else
        conn(i+1,2) = 0;
    end
end
figure
semilogy(conn(:,1),conn(:,2));
hold on
a=semilogy(conn(:,1),conn(:,2),'k.');
set(a,'MarkerSize',10);
xlim([0 1]);
ylim([0.00001 1]);
xlabel("cut-off");
ylabel("connectivity");
strValues = strtrim(cellstr(num2str(conn(:,2))));
text(conn(:,1),conn(:,2),strValues,'VerticalAlignment','bottom');

hold off
drawnow
savefig(output + "_connectivity.fig")
writematrix (conn, output + "_connectivity.txt");

allcorr = size(C(C>threshcorr))/2;

main_analysis(28,1) = "all correlations above threshcorr";
main_analysis(28,2) = allcorr(1);
main_analysis(29,1) = "connectivity (all correlations above threshcorr / corrleft)";
main_analysis(29,2) = allcorr(1) / corrleft;
end