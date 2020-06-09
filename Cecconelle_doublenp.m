%% _____double np another way
clc; close all; clear;
tic;
 % Neutron energy
 Eneutron=2.45; % MeV
 
 

 % Define parameters Beta
% clear('prompt');
 % prompt= 'Input Beta. Cecconell: 1.5 |Verbinski: 1.23 | EJ299: 1.63 (1.44| 1.82) . >>>';
 % Beta=str2num(input(prompt, 's' ));
Beta= 1.5;
 % Define parameter k
 %clear('prompt');
 % prompt= 'Input k. Cecconell: 0.21 |Verbinski: 0.265 |EJ299: 0.12 (0.08 | 0.16). >>>';
 %k=str2num(input(prompt, 's' ));
k=0.21;
%     k=0.205; % Cecconello: MeV^(-Beta) k=0.204; k ~ 0.21MeV^(-Beta)

% Define maximum light output 
 MaxLO=k*(Eneutron)^(Beta);

% let the oberserved light output be L
%     L=MaxLO;
    L=round(linspace(0, MaxLO, 1e5), 5);
    %iL=1:1:length(L);

% find the index correspording the observed light output and set the array
    % Define parameter k
%     clear('prompt');
%     prompt= 'What is the observed Light output? 0.56945 >>>';
%     L_Obs=str2num(input(prompt, 's' ));

    L_Obs=0.56945;
    
    iL=find(L==L_Obs);
    clear('L_Obs');
    
    % ___L1___
    clear('Ep1_dnp', 'LO1_dnp')
    
    % Ep___
%     Ep1_dnp= linspace(0, (L(iL)/k)^(1/Beta), 1e3);  % Ep1
%     LO1_dnp=k.*(Ep1_dnp).^(Beta);  % LO1
    % Ep___
    
    LO1_dnp=linspace(0, L(iL), 1e7);
    Ep1_dnp= (LO1_dnp./k).^(1/Beta);

    % ___L2___
    clear('Ep2_dnp', 'LO2_dnp')
%     LO2_dnp=L(iL) - LO1_dnp;   %LO2
%     Ep2_dnp= (LO2_dnp ./ k).^(1/Beta);  %Ep2
    
    % Ep___
    LO2_dnp=L(iL)- LO1_dnp;   %LO2
    Ep2_dnp= (LO2_dnp ./ k).^(1/Beta);  %Ep2
    % Ep___
    
    
     clear('Ed_dnp')
     Ed_dnp=(Ep2_dnp + Ep1_dnp); 
    

   
    % PDF
    PDF_L1_dnp= (1/Eneutron).*(1./(Beta.*LO1_dnp)).*((LO1_dnp./k).^(1/Beta));
    
    PDF_L2_dnp =(1./(Eneutron-Ep1_dnp)).*(1./(Beta.*LO2_dnp)).*((LO2_dnp./k).^(1/Beta));
    
    
    
    %  clear('PDF_resultant_dnp')
    PDF_resultant_dnp = (PDF_L1_dnp.*PDF_L2_dnp);
    
    


     
     clear("i_EdLTEn");
     
     i_EdLTEn=find((Ed_dnp-Eneutron)>0);
     if length(i_EdLTEn)>0
         PDF_resultant_dnp(i_EdLTEn)=0;
         PDF_L1_dnp(i_EdLTEn)=0;
         PDF_L2_dnp(i_EdLTEn)=0;
     end
     
     clear("i_EdLTEn");

     
     % set Nan ==0
      PDF_resultant_dnp(isnan(PDF_resultant_dnp))=0;
     PDF_L1_dnp(isnan(PDF_L1_dnp))=0;
     PDF_L2_dnp(isnan(PDF_L1_dnp))=0;
   
% convolution of PDF_L1_dnp with PDF_L2_dnp
        clear('U', 'V', 'W', 'm','n');
        x=fliplr(PDF_L1_dnp); %_1
        h=fliplr(PDF_L2_dnp); %_1
        
        
%           x=PDF_L1_c;            %_2
%           h=fliplr(PDF_L2_c);    %_2

%           x=PDF_L1_c;    %_3
%           h=PDF_L2_c;    %_3

%           h=PDF_L1_c;    %_4
%           x=PDF_L2_c;    %_4
%           
%           h=PDF_L1_c;            %_5
%           x=fliplr(PDF_L2_c);    %_5

%         x=PDF_Ep_1;
%         h=PDF_Ep_2;

        m=length(x);
        n=length(h);
       
        U=[x,zeros(1,n)]; 
        V=[h,zeros(1,m)]; 

		 % Make sure all Nan component ==0
        if length(find(isnan(U)==1))>0
            U(find(isnan(U)==1))=0;
        end
        
        if length(find(isnan(V)==1))>0
            V(find(isnan(V)==1))=0;
        end

            nR=1;
			nC=1;
			
%             
        for nR=1:1:m+n-1
            

           
            nC=1;
			totL(nR)=0;
            for nC=1:1:nR
             %   W(nR,nC)= U(nC)*V(nR-nC+1);
			totL(nR)=totL(nR)+U(nC)*V(nR-nC+1);

                nC=nC+1;
            end

            nR=nR+1;
            
            if mod(nR, 10000) ==0
            sprintf('No of row= %s / %s', string(nR), string(length(V)))
            end

        end
        
        clear('nC','nR','U','V')
           % Ep__ 
%         Epdnp_x= Ep1_dnp(1):Ep1_dnp(1):Ep1_dnp(1)*length(dnp_nY);
%         LOdnp_x= k*(Epdnp_x).^(Beta);
        % Ep__
        
        %LO__
        LO_x_step=LO1_dnp(2)-LO1_dnp(1);
        LOdnp_x= LO_x_step:LO_x_step:LO_x_step*length(totL);
        %LO__
        
        %___ To modified the  W(nR,nC)= U(nC)*V(nR-nC+1); _____
%          W(isnan(W)==1)=0;
%         dnp_nY=sum(W, 2);
%         clear('i')
%         i=find(LOdnp_x> round(k*(Eneutron)^(Beta), 1));
%         dnp_nY(i)=0;
%         clear('i')
        
        %___ To modified the totL(nR)=totL(nR)+U(nC)*V(nR-nC+1);
        clear('i')
        i=find(LOdnp_x> round(k*(Eneutron)^(Beta), 1));
        totL(i)=0;
        clear('i');
      
        % plot the Ed=Ep1 + Ep2; P(L2|L1)*p(L1) and convolution data
        figure%('Visible', 'off');
        h1=subplot(1,3,1)
        plot(LO1_dnp, Ed_dnp, 'g:', 'LineWidth',3); hold on;
        plot(LO1_dnp, Ep1_dnp, 'r--', 'LineWidth',3); hold on;
        plot(LO1_dnp, Ep2_dnp, 'b--', 'LineWidth',3)
        line([0, 0.8], [Eneutron Eneutron]); hold on;
%         ylim([0 3])
        xlabel('LO1')
        ylabel('Total energy (MeV)')
        legend('Ed = Ep1 + Ep2','Ep1','Ep2', 'Eneutron= 2.45 MeV', 'Location', 'southoutside');
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2);
        pbaspect([1.5 1 1]);
        axis square;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        
        h2=subplot(1,3,2)
        area(LO1_dnp,PDF_resultant_dnp, 'FaceColor', [0.6 0.6 0.6] ); hold on;
        plot(LO1_dnp, PDF_L1_dnp,'b-', 'LineWidth', 1.5 ); hold on;
        plot(LO1_dnp, PDF_L2_dnp,'r-', 'LineWidth', 1.5  ); hold on;
        str=strcat('Observed L= ', num2str(L(iL)), ' Area= ', num2str(trapz(LO1_dnp, PDF_resultant_dnp)));
        xlabel('LO of 1st recoil proton (MeVee)')
        ylabel('PDF_p(L)(a.u.)')
        axis square;
        title(str)
        xlim([0 L(iL)+0.1]); 
    %     ylim([0 13]);
        grid on; 
        grid minor;
    %     ylim([0 0.6]);
        legend('p(L)','p(L1)','p(L2|L1)', 'Location', 'southoutside')
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2);
        pbaspect([1.5 1 1]);
        axis square
        
        h3=subplot(1,3,3)
        plot(LOdnp_x, totL./trapz(LOdnp_x, totL), 'k-', 'LineWidth', 2); hold on;
%         plot(LOdnp_x, dnp_nY./trapz(LOdnp_x, dnp_nY'), 'k-', 'LineWidth', 2);
        xlabel('LO');
        ylabel('PDF P(L2|L1)P(L1)');
        xlim([0 1.0]);
        ylim([0 3.5]);
        clear('str');
        str={strcat('Observed L = ', num2str(L(iL)), '| k = ', num2str(k), '| Beta = ', num2str(Beta));
              strcat('steps= ', num2str(length(LO1_dnp)))};
        title(str)
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2);
        legend('p(L)=p(L2|L1){\otimes}p(L1)', 'Location', 'southoutside')
        pbaspect([1.5 1 1]);
        axis square;
		
		  clear('fN');
        fN=sprintf('dnp_L_%s_k_%s_Beta_%s_Steps_%s.png', num2str(L(iL)), num2str(k), num2str(Beta), num2str(length(LO1_dnp)));
        saveas(gcf, fN, 'png'); clear('fN');
       
        
        figure;
%         plot(LOdnp_x, dnp_nY./trapz(LOdnp_x, dnp_nY'), 'k-', 'LineWidth', 2);hold on;
%         plot(LOdnp_x, totL./trapz(LOdnp_x, totL), 'r:', 'LineWidth', 2);hold on;
%         plot(LOdnp_x, dnp_nY./trapz(LOdnp_x, dnp_nY'), 'k-', 'LineWidth', 2); hold on;
        plot(LOdnp_x, totL./trapz(LOdnp_x, totL), 'k-', 'LineWidth', 2); hold on;
        xlabel('LO');
        ylabel('PDF P(L2|L1)P(L1)');
        xlim([0 1.0]);
%         ylim([0 3.5]);
        clear('str');
        str={strcat('Observed L = ', num2str(L(iL)), '| k = ', num2str(k), '| Beta = ', num2str(Beta));
            strcat('steps= ', num2str(length(LO1_dnp)))};
        title(str)
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2);
        legend('p(L)=p(L2|L1){\otimes}p(L1)', 'Location', 'southoutside')
        pbaspect([1.5 1 1]);
        
        clear('fN');
        fN=sprintf('dnp_convL_%s_k_%s_Beta_%s_Steps_%s.png', num2str(L(iL)), num2str(k), num2str(Beta), num2str(length(LO1_dnp)));
        saveas(gcf, fN, 'png'); clear('fN');
        

        % output convolution result
        LOdnp_x=LOdnp_x';
        totL=totL';
        tbl=table(round(LOdnp_x, 5), round(totL, 5));
        clear('fN');
        fN=sprintf('data_convL_%s_k_%s_Beta_%s_Steps_%s.txt', num2str(L(iL)), num2str(k), num2str(Beta), num2str(length(LO1_dnp)));
        writetable(tbl, fN)
%         close all;
        
      
%     iL=iL+1;
%     end
    toc;