%load gaimc software package from matlab
%(https://de.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code?s_tid=srchtitle)
%David Gleich (2022). gaimc : Graph Algorithms In Matlab Code (https://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code), MATLAB Central File Exchange.

function [main_analysis, CC] = networktheory(C, main_analysis, threshcorr)

%simplifying data and creating random network data with same network size
%and number of connections thus preserving the same degree on average
CC = zeros(size(C,1));
x = zeros(size(C,1));
for i=1 : size(C,1)
    for ii=1 : size(C,1)
        if C(ii,i) > threshcorr && x(i,ii) == 0
            CC(ii,i) = 1;
            CC(i,ii) = 1; %Matrix muss symmetrisch sein!
        end
        x(ii,i) = 1;
    end
end

conncells=0;
for z=1 : size(CC,2)
    if sum(CC(:,z)) > 0
        conncells=conncells+1; %Zählt Zellen die mindestens eine connection haben
    end
end

randCC = zeros(size(CC,2)); 

for i=1 : sum(sum(CC))/2 %es wird die gleiche Anzahl von peaks wie in der original matrix verteilt
    y=0;
    while y==0
        randcell1 = 1+floor(conncells * rand); %damit in der random matrix nur so viele Zellen verbunden sind wie auch in der original matrix auch Verbunden sind
        randcell2 = 1+floor(conncells * rand);
        if randCC(randcell1,randcell2) == 0 && randcell1 ~= randcell2 %while schleife wird so lange wiederholt bis eine Kombination gefunden wurde, die noch nicht benutzt wurde
            y=1;
        end
    end
    randCC(randcell1,randcell2) = 1;
    randCC(randcell2,randcell1) = 1;
end

%mean shortest path (SP)
CCgraph = graph(CC);
SP = distances(CCgraph);
SP(SP==Inf) = NaN;
meanSP = mean(SP,'all','omitnan');

randCCgraph = graph(randCC);
randSP = distances(randCCgraph);
randSP(randSP==Inf) = NaN;
randmeanSP = mean(randSP,'all','omitnan');

main_analysis(30,1) = "mean shortest path (original network)";
main_analysis(30,2) = meanSP;
main_analysis(31,1) = "mean shortest path (random network)";
main_analysis(31,2) = randmeanSP;


%mean clustering coeficcient (Cl)
sparseCC = sparse(CC);
Cl = clustercoeffs(sparseCC);
meanCL = mean(Cl);

randsparseCC = sparse(randCC);
randCl = clustercoeffs(randsparseCC);
randmeanCL = mean(randCl);

main_analysis(32,1) = "mean clustering coeficcient (original network)";
main_analysis(32,2) = meanCL;
main_analysis(33,1) = "mean clustering coeficcient (random network)";
main_analysis(33,2) = randmeanCL;


%calculate the degree (they will always have same value)
%degree = (sum(sum(CC)))/size(C,1)
%randdegree = (sum(sum(randCC)))/size(C,1)

end

