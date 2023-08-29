function numberofconnections(C, thresh, output, periodicity)
Cnum = zeros(size(C,1),1);
Cnum_perio = zeros(size(C,1),1);
not = zeros(size(C,1));

for i=1 : size(C,1)
    for ii=1 : size(C,1)
        if C(i,ii) >= thresh && not(i,ii) == 0
            not(ii,i) = 1;
            Cnum(i,1)=Cnum(i,1)+1;
            Cnum(ii,1)=Cnum(ii,1)+1;
            if periodicity(i) == 1
                Cnum_perio(i,1)=Cnum_perio(i,1)+1;
            end
            if periodicity(ii) == 1
                Cnum_perio(ii,1)=Cnum_perio(ii,1)+1;
            end
        end
    end
end
%proablility distribution P(k):
propdist = zeros(max(Cnum),3);
for j=1 : max(Cnum)
    m=size(Cnum(Cnum==j));
    m_perio=size(Cnum_perio(Cnum_perio==j));
    propdist(j,1)= j;
    propdist(j,2) = m(1);
    propdist(j,3) = m_perio(1);
end
writematrix (propdist, output + "_propablility distribution.txt");

figure
histogram(Cnum);
grid on
xlabel("number of connections/cell");
drawnow
savefig(output + "_propablility distribution.fig")
end