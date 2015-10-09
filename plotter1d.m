
% plotter1d.m
% By: Devin Light
% ------

clc;
clear all;

set(0,'defaultfigureposition',[180 520 3*180 2*180],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',18); format compact, format long
%%
cd('/Users/Devin/Desktop/R/NodalDG/nodal_final');

tests = { ...
         'sqwave', ... % square wave
         'cos', ... % cos
         'cos2', ... % cos**2
         'cos4', ... % cos**4
         '2wave', ... % 2+sin(6*pi*x)+sin(8*pi*x)
         };
methods = { ...
          'nod', ...
          'zshu',...
          };
      
res = {'1','2','3','4','5'};
which_res = res(3);     
which_test = tests(1);
ncfilename = strcat('dg1d_' ,which_test{1},'.nc');

%% Cycle through methods

stat        = 2;
which_res   = res(3);
which_test  = tests(1);
printExtrema = 1;

for n=1:4
    if(n==1) 
        methname = 'Unlimited';
        nc = ['_ndgunlim/' ncfilename];
        file = ['figures/nodal/nod_' which_test{1}];
        nodal_unlim = plot_1dadv(methname,nc,which_res,file,stat,printExtrema);
    elseif(n==2)
        methname = 'ZS Limited';
        nc = ['_ndgzhshu/' ncfilename];
        file = ['figures/zshu/zshu_' which_test{1}];
        nodal_shu = plot_1dadv(methname,nc,which_res,file,stat,printExtrema);
    elseif(n==3)
        methname = 'TMAR Limited';
        nc =['_matrunc/' ncfilename];
        nodal_ma = plot_1dadv(methname,nc,which_res,file,stat,printExtrema);
    elseif(n==4)
        methname = 'TMAR FCT';
        nc =['_tmarFCT/' ncfilename];
        nodal_fct = plot_1dadv(methname,nc,which_res,file,stat,printExtrema);
    end
end

%% Animate movie
f = figure();
hold = nodal_shu;
tmp = squeeze(hold.data(1,:));
ics = squeeze(hold.data(1,:));
fin = squeeze(hold.data(end,:));
axis([0 1 -.1 1.2]);

for i=1:length(hold.t)
    tmp = squeeze(hold.data(i,:));
    ics = squeeze(hold.data(1,:));

    loc = tmp<0.0;
    h = plot(hold.x,tmp,'k.-',hold.x,ics,'r.-',hold.x(loc),tmp(loc),'gx');
    axis([0 1 -.1 1.1]);
    ftitle = ['t=',num2str(hold.t(i))];
    title(ftitle);

    pause(0.2);
end
%}
%% Make revolution comparison figures
f = figure();
hold = nodal_shu; x = hold.x;
ics = squeeze(hold.data(1,:));

t1rev = hold.t == 1.0;
rev1 = squeeze(hold.data(t1rev,:));

t10rev = hold.t == 10.0;
rev10 = squeeze(hold.data(t10rev,:));

loc = rev10 < 0.0;
plot(x,ics,'r.-',x,rev1,'k.-',x,rev10,'b.-',x(loc),rev10(loc),'gx');
axis([0 1 -.1 1.2]);
leg1 = 'ICs';
leg2 = ['t=',num2str(hold.t(t1rev))];
leg3 = ['t=',num2str(hold.t(t10rev))];
legend(leg1,leg2,leg3);
ftitle = [hold.method ' Advection at 3 times on [' num2str(round(x(1))) ',' num2str(round(x(end))) '] domain'];
title(ftitle);


name = strcat('ndgshu_longadv.pdf');
print(f,'-dpdf',name)
%%
figpos = [180 520 2*180 2*180];

hold = Q2;
xHold = x2;

ics = squeeze(hold(1,:));
final = squeeze(hold(end,:));

fig = figure('Position',figpos);
subaxis (1,1,1, 'Padding', 0, 'MarginRight', 0.05,'MarginLeft', 0.05,'MarginTop', 0.05,'MarginBottom', 0.075);
plot(xHold,ics,'r-',xHold,final,'k.-');
axis([0 1 -.1 1.2]);

set(fig,'units','centimeters');
op = get(fig,'OuterPosition');
set(fig, 'PaperUnits','centimeters'); set(fig, 'PaperSize', [op(3) op(4)]);
set(fig, 'PaperPositionMode', 'manual'); set(fig, 'PaperPosition',[0 0  op(3) op(4)]);
saveas(fig, 'sqwave_1d_unmod.pdf', 'pdf');

hold = Q3;
xHold = x3;

ics = squeeze(hold(1,:));
final = squeeze(hold(end,:));

fig = figure('Position',figpos);
subaxis (1,1,1, 'Padding', 0, 'MarginRight', 0.05,'MarginLeft', 0.05,'MarginTop', 0.05,'MarginBottom', 0.075);
plot(xHold,ics,'r-',xHold,final,'k.-');
axis([0 1 -.1 1.2]);

set(fig,'units','centimeters');
op = get(fig,'OuterPosition');
set(fig, 'PaperUnits','centimeters'); set(fig, 'PaperSize', [op(3) op(4)]);
set(fig, 'PaperPositionMode', 'manual'); set(fig, 'PaperPosition',[0 0  op(3) op(4)]);
saveas(fig, 'sqwave_1d_PDDG.pdf', 'pdf');

hold = Q4;
xHold = x4;

ics = squeeze(hold(1,:));
final = squeeze(hold(end,:));

fig = figure('Position',figpos);
subaxis (1,1,1, 'Padding', 0, 'MarginRight', 0.05,'MarginLeft', 0.05,'MarginTop', 0.05,'MarginBottom', 0.075);
plot(xHold,ics,'r-',xHold,final,'k.-');
axis([0 1 -.1 1.2]);

set(fig,'units','centimeters');
op = get(fig,'OuterPosition');
set(fig, 'PaperUnits','centimeters'); set(fig, 'PaperSize', [op(3) op(4)]);
set(fig, 'PaperPositionMode', 'manual'); set(fig, 'PaperPosition',[0 0  op(3) op(4)]);
saveas(fig, 'sqwave_1d_zshu.pdf', 'pdf');


%% Make plot with interpolating polynomial on top of nodal values
nodalVals = squeeze(nodal_shu.data(2,:));
nodalICs = squeeze(nodal_shu.data(1,:));
quadWeights = nodal_shu.weights;
nDegree = nodal_shu.N; xQuad = nodal_shu.x;

nelem = length(nodalVals)/(nDegree+1);
xPlot = linspace(0,1,length(nodalVals)); qPlot = []; massJ = [];
for j=1:nelem
    elemVals = nodalVals(1+(j-1)*(nDegree+1):j*(nDegree+1));
    elemNodes = xQuad(1+(j-1)*(nDegree+1):j*(nDegree+1));
    p = polyfit(elemNodes',elemVals,nDegree);
    qPlot = [qPlot polyval(p,xPlot(1+(j-1)*(nDegree+1):j*(nDegree+1)))];
    
    quadVals = nodalVals(1+(j-1)*(nDegree+1):j*(nDegree+1));
    massJ = [massJ (0.5/nelem)*sum(quadWeights.*quadVals')];
end
totalMass = sum(massJ);
exMass = 0.2;
massError = abs(exMass-totalMass);

exact = squeeze(nodal_shu.data(1,:));
plot(xPlot,qPlot,'g.-',xQuad,exact,'r--',xQuad,nodalVals,'ko-')
