clearvars;

pkloc = '/Users/davidhausmann/Desktop/Neuer Ordner/190314_100000_S24 wt_rhod2_72h after seeding_calcium time series_well2b_1_Krebs solution_pklocs.txt';
pfile = readmatrix(pkloc);

[frame,cellnumber] = find(pfile);

%dotplot from pkloc
figure
b = plot(frame,-cellnumber,'k.');
set(b,'MarkerSize',5);
xlabel("Frame");
ylabel("Cellnumber");
xlim([0 max(frame)]);
ylim([-max(cellnumber) 0]);

%densityplot from pkloc
figure
h=histogram2(frame,-cellnumber,'DisplayStyle','tile');
h.BinWidth = [35 15];
xlabel("Frames");
ylabel("# of cell");