FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(poutsf2(:,1)),poutsf2(:,9),'m--^',log10(poutsf2(:,1)),poutsf2(:,6),'b--o',log10(poutm(:,1)),poutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless pore pressure P/P_0'); 
legend('NS-FEM','FEM-T3',' Manoharan (1995)')
title('Pore pressure with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');


FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(uoutsf2(:,1)),uoutsf2(:,9),'m--^',log10(uoutsf2(:,1)),uoutsf2(:,6),'b--o',log10(uoutm(:,1)),uoutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless settlement 100u/a'); 
legend('NS-FEM','FEM-T3','Manoharan (1995)')
title('Settlement with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');




f_out1=3*200*linspace(1,30,30)/30;
figure
plot(abs(u_out),abs(f_out1),'--b*','linewidth',2);
xlabel({'Displacement (m)'},'FontSize',16);
ylabel({'Vertical force (KN)'},'FontSize',16);
title('F-U curve (load control)')



FontSize=18;
figure1=figure(1);
u_out=e2pfooting;
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(abs(u_out(:,1)),abs(f_out1),'m--^',abs(u_out(:,2)),abs(f_out1),'b--o',abs(u_out(:,3)),abs(f_out1),'r--*')
xlabel('Displacement (m)');
ylabel('Vertical force (KN)'); 
legend('FEM-T3','NS-FEM','NS-PFEM')
title('F-U curve (load control)');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');


nsteps=20;
step_array=linspace(1,nsteps,nsteps);
FontSize=18;
figure1=figure(1);
u_out=-e2pslope;
f_out1=step_array;
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(abs(u_out(:,1)),abs(f_out1),'m--^',abs(u_out(:,2)),abs(f_out1),'b--o',abs(u_out(:,3)),abs(f_out1),'r--*')
xlabel('Displacement (m)');
ylabel('Load step'); 
legend('FEM-T3','NS-FEM','NS-PFEM')
title(' Top displacement (slope)');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');

t1=creep5(:,1);
t2=creep5(:,2);
t3=creep5(:,3);

FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(t1),-t3,'b-*');
xlabel('Time lg(T (h))');
ylabel('Pore pressure P(kPa)'); 
title('Pore pressure with time');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');

FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(t1),t2,'r-*');
xlabel('Time lg(T (h))');
ylabel('Settlement(m)'); 
title('Settlement with time');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');












FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(SNSpout(:,1)),SNSpout(:,3),'m--^',log10(SNSpout(:,1)),SNSpout(:,2),'b--o',log10(SNSpout(:,1)),SNSpout(:,4),'k--d',log10(poutm(:,1)),poutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless pore pressure P/P_0'); 
legend('IS+PS','IS','PS',' Manoharan (1995)')
title('Pore pressure with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');


FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(SNSuout(:,1)),SNSuout(:,3),'m--^',log10(SNSuout(:,1)),SNSuout(:,2),'b--o',log10(SNSuout(:,1)),SNSuout(:,4),'k--d',log10(uoutm(:,1)),uoutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless settlement 100u/a'); 
legend('IS+PS','IS','PS',' Manoharan (1995)')
title('Settlement with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');


FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(SNSpout(:,1)),SNSpout(:,3),'b--o',log10(poutm(:,1)),poutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless pore pressure P/P_0'); 
legend( 'SNS-FEM',' Manoharan (1995)')
title('Pore pressure with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');
box off


FontSize=18;
figure1=figure(1);
set(1,'color','w');
set(1,'Position',[500 200 600 500]);
plot(log10(SNSuout(:,1)),SNSuout(:,3)/1.01,'b--o',log10(uoutm(:,1)),uoutm(:,2),'r--*')
xlabel('Dimensionless time lg(T_v)');
ylabel('Dimensionless settlement 100u/a'); 
legend('SNS-FEM',' Manoharan (1995)')
title('Settlement with T_v');
set(gca,'Fontsize',FontSize,'FontName','Times new Roman');
set(get(gca,'XLabel'),'FontSize',FontSize,'FontName','Times new Roman');
set(get(gca,'YLabel'),'FontSize',FontSize,'FontName','Times new Roman');
box off