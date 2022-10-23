function main_analysis = networkgraphs(hubs, periodicity, CC, output, do_periodicity, main_analysis)
%plots the original network and the original network without periodic cells
%and calculates the larges cluster

node = 4;
edge = 1.65;

hh = NaN(1,size(hubs,2));
for z=1 : size(hubs,2)
    if hubs(z) == 1
        hh(z) = z;
    end
end

pp = NaN(1,size(periodicity,2));
for z=1 : size(periodicity,2)
    if periodicity(z) == 1
        pp(z) = z;
    end
end

stop = size(CC,1);
CC_2=CC;
kk=0;
for k=1 : stop %löscht alle Zellen, die nicht zum Netzwerk gehören
    kk=kk+1;
    if sum(CC(:,k),'omitnan') == 0
        CC_2(kk,:)=[];
        CC_2(:,kk)=[];
        hh(1,kk:end)=hh(1,kk:end)-1;
        hh(:,kk)=[];
        pp(1,kk:end)=pp(1,kk:end)-1;
        pp(:,kk)=[];
        kk=kk-1;
    end
end

pp=pp(~isnan(pp));
pp_size = size(pp(1,:));

hh=hh(~isnan(hh));
CC_2 = sparse(CC_2);

figure
G=graph(CC_2);
g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);   %g=plot(G,'NodeColor','k');
layout(g, 'force','UseGravity',true)
highlight(g, hh(1,:),'NodeColor','#AA0000');
drawnow
output_figure = sprintf('%s_network_hubs.fig', output);
savefig(output_figure);
exportgraphics(gcf,output + "_network_hubs.pdf",'ContentType','vector')


figure
G=graph(CC_2);
g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
layout(g, 'force','UseGravity',true)
highlight(g, pp(1,:),'NodeColor','#AA0000');
drawnow
output_figure = sprintf('%s_network_perio.fig', output);
savefig(output_figure);
exportgraphics(gcf,output + "_network_perio.pdf",'ContentType','vector')


figure
G=graph(CC_2);
g=plot(G, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge, 'NodeColor', [0 0.447 0.741]);
layout(g, 'force','UseGravity',true)
highlight(g, pp(1,:),'NodeColor','#00ff00', 'MarkerSize', node*2.5);
hold
%g2=plot(G, 'NodeLabel',{}, 'Marker','none', 'Edgecolor','none');
g2=plot(G, 'NodeLabel',{}, 'NodeColor', [0 0.447 0.741],'Edgecolor','none', 'MarkerSize', node);
layout(g2, 'force','UseGravity',true)
highlight(g2, hh(1,:),'NodeColor','#AA0000');
drawnow
output_figure = sprintf('%s_network_hubsandperio.fig', output);
savefig(output_figure);
exportgraphics(gcf,output + "_network_hubsandperio.pdf",'ContentType','vector')


if do_periodicity > 0
    periodicity2=periodicity;
    CC_notperio = full(CC_2);
    for k=1 : size(CC_notperio,1)
        for l=1 : size(CC_notperio,1)
            if ismember(k,pp) == 1
                periodicity2(k) = 1;
            else
                periodicity2(k) = 0;
            end
            if ismember(l,pp) == 1
                periodicity2(l) = 1;
            else
                periodicity2(l) = 0;
            end
            if periodicity2(k) == 1 || periodicity2(l) == 1
                CC_notperio(k,l) = 0;
                CC_notperio(l,k) = 0;
            end
        end
    end
    
    CC_notperio = sparse(CC_notperio);
    figure
    H=graph(CC_notperio);
    h=plot(H, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
    layout(h, 'force','UseGravity',true)
    highlight(h, pp(1,:),'NodeColor','#646464');
    
    
    drawnow
    output_figure = sprintf('%s_network_notperio.fig', output);
    savefig(output_figure);
    exportgraphics(gcf,output + "_network_notperio.pdf",'ContentType','vector')

    
    
   
    rando=periodicity;
    CC_rando = full(CC_2);
    
    if size(CC_rando,2) >= sum(periodicity,'all','omitnan')
        r = randperm(size(CC_rando,2),pp_size(2));
        for k=1 : size(CC_2,1)
            for l=1 : size(CC_2,1)
                if ismember(k,r) == 1
                    rando(k) = 1;
                else
                    rando(k) = 0;
                end
                if ismember(l,r) == 1
                    rando(l) = 1;
                else
                    rando(l) = 0;
                end
                if rando(k) == 1 || rando(l) == 1
                    CC_rando(k,l) = 0;
                    CC_rando(l,k) = 0;
                end
            end
        end

        [~,binsizes_rando] = conncomp(graph(CC_rando));
        largestcluster_rand = max(binsizes_rando);
        
        CC_rando = sparse(CC_rando);
        figure
        H=graph(CC_rando);
        h=plot(H, 'NodeLabel',{}, 'MarkerSize', node, 'LineWidth', edge);
        layout(h, 'force','UseGravity',true);
        highlight(h, r,'NodeColor','#646464');

        drawnow
        
        output_figure = sprintf('%s_network_notrand_%.0f.fig', output,largestcluster_rand);
        savefig(output_figure);
        output_figure = sprintf('%s_network_notrand_%.0f.pdf', output,largestcluster_rand);
        exportgraphics(gcf,output_figure,'ContentType','vector')

    end
    
    [~,binsizes] = conncomp(graph(CC_2));
    largestcluster = max(binsizes);
    [~,binsizes_notperio] = conncomp(graph(CC_notperio));
    largestcluster_perio = max(binsizes_notperio);
    [~,binsizes_rando] = conncomp(graph(CC_rando));
    largestcluster_rand = max(binsizes_rando);
    
    main_analysis(34,1) = "Size of largest network";
    if sum(largestcluster,'all','omitnan') > 0
        main_analysis(34,2) = largestcluster;
    else
        main_analysis(34,2) = 0;
    end
    
    main_analysis(35,1) = "Size of largest network after removing all periodic cells";
    if sum(largestcluster_perio,'all','omitnan') > 0
        main_analysis(35,2) = largestcluster_perio;
    else
        main_analysis(35,2) = 0;
    end
    
    main_analysis(36,1) = "Size of largest network after removing random cells in the same number as periodic cells";
    if sum(largestcluster_rand,'all','omitnan') > 0
        main_analysis(36,2) = largestcluster_rand;
    else
        main_analysis(36,2) = 0;
    end
    
    
    
    
    
end