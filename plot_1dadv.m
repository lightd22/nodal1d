% plot_1dadv.m
% Reads in netCDF files and plots data for 1d advection
% By: Devin Light
% ---

function [out] = plot_1dadv(methname,nc,res,file,stat,print_extrema)
    
    Qname = strcat('q',res{1});
    xname = strcat('x',res{1});
    muname = strcat('mu',res{1});
    weightsname = 'qweights';
    nodesname = 'qnodes';
    tname = 'time';
    
    out.data = nc_varget(nc, Qname);
    out.x = nc_varget(nc, xname);
    out.t = nc_varget(nc, tname);
    out.weights = nc_varget(nc,weightsname);
    out.nodes = nc_varget(nc,nodesname);
    out.mu = nc_varget(nc,muname);
    out.N = length(out.nodes) - 1;
    out.method = methname;
    
    if(stat == 1)
    elseif(stat == 0)
    else
        fig = figure();
        h = subplot(1,1,1);
        tmp = squeeze(out.data(end,:));
        tmp_ics = squeeze(out.data(1,:));
        x = out.x; t = out.t;
        plot(x,tmp,'k.-',x,tmp_ics,'r-.');
        refline(0,0)

if(print_extrema == 1)
        text(0.65, .75,'Max','FontSize',14);
        text(0.65, .68,['Exact:',num2str(max(tmp_ics))],'FontSize',12);
        text(0.65, .58,['Apprx:',num2str(max(tmp))],'FontSize',12);

        text(0.20, .75,'Min','FontSize',14);
        text(0.20, .68,['Exact:',num2str(min(tmp_ics))],'FontSize',12);
        text(0.20, .58,['Apprx:',num2str(min(tmp))],'FontSize',12);
end
        
        ftitle = {strcat(methname, sprintf(', Time: %0.2f sec ',t(end))),[' nelem=' num2str(size(tmp,2)/(out.N + 1))]};
        title(ftitle);
        axis square;
        xlim([0 1]);
        ylim([-0.1 1.1]);

        name = strcat(file,'.pdf');
        print(fig,'-dpdf',name)

    end
end