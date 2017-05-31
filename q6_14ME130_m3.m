%--------------------------------------------------------------------------------------------------%
%--MATLAB-Individual Assignment-MA03 {Individual Question-Assignment}--%
%--Q6-Slotted Link QR Mechanism-Position,Vector,Acceleration,Static,Dynamic,Shaking Force Analysis%
%--Isac Rajan--%
%--14ME130--%
%--Group M--%
%--------------------------------------------------------------------------------------------------%

%Input dialog box for taking input from User%
msg=msgbox('Please NOTE: DFA, SFA stands for Dynamic Force Analysis and Shaking Force Analysis respectively.');
choice=menu('Default values','Set 1 + DFA + SFA','Set 2 + DFA + SFA','Set 3 + DFA + SFA','Snapshot + Static Force Analysis');
prompt={'Length of frame', 'Length of input crank', 'Length of slotted lever',...
    'Length of output link', 'Height of slider from ground','Initial Angular Velocity'...
    ,'Angular acceleration','Force at Slider(newtons)','mass of slider','Linear mass density(kg/units)'};
title='Input for Crank Slotted Mechanism';
lines=[1 60];
if choice==1
    default={'5', '3', '10', '5', '6', '1', '0.05','10','2','0.1'};
elseif choice==2
    default={'15', '9', '30', '15', '18', '1.2', '0.06','30','20','0.1'};
elseif choice==3
    default={'500', '300', '1000', '500', '600', '0.8', '0.07','1250','150','0.1'};
elseif choice==4 %if choice 4 is chosen Snapshot is given%
    ss=inputdlg({'Length of frame', 'Length of input crank', 'Length of slotted lever',...
        'Length of output link', 'Height of slider from ground','Crank-Angle(Degrees)',...
        'Angular Velocity','Angular acceleration', 'Force - F (NEWTONS)'},'Snapshot+SFA',[1 40],{'5', '3', '10',...
        '5', '6','60','1','1','10'},'on');
end
if (choice~=0 && choice~=4)
    dlg=inputdlg(prompt,title,lines,default,'on');
%Converting the string input to double precision numerals
l1=str2double(dlg(1));
l2=str2double(dlg(2));
l4=str2double(dlg(3));
l5=str2double(dlg(4));
l7=str2double(dlg(5));
w0=str2double(dlg(6));
aa=str2double(dlg(7));
F=str2double(dlg(8));
m=str2double(dlg(9));
rho=str2double(dlg(10));
%Checking validity of link lengths
if (l4 < 2*l2) || (((l4 - l5) > l7) || (l5 > l4)) ||...
        (l1 + l2 > l4) || (l7 > 0.707*l4 + l5) || (l2>l1) || l2<=0 || l1<=0 || l4<=0 || l5<=0 || l7<=0
    h = errordlg('The links are not valid to form a mechanism. The mechanism cannot be simulated','Error');
    return
else
    w=0; t2=0; t3=0; %initialization
    h=figure('Name','Crank Slotted Lever Mechanism - Simulation & Analysis');
    t=5; %a random value for time to start with
    gg=msgbox('Please see the required TORQUE and SHAKING FORCE values on the COMMAND PROMPT');
    while ishandle(h)
        %Position Analysis and Configuration Diagram of Crank Slotted Lever Mechanism%
        subplot(2,2,1);
        text(l1/10,0, 'O2');
        plot([0 0],[0 l1], '--o','linewidth',1.5);
        hold off; 
        text(-l1/3,l1, 'O4');
        axis equal;
        axis([-(3*l5) (3*l5) -l1/2 (l1+l4)]);
        pause(0.05);
        t2=(w0*t + 0.5*aa*(t*t));
        t3=pi+atan( (l2*sin(t2)+l1)/(l2*cos(t2)) );
        l3=l2 * cos(t2) / cos(t3);
        plot([0 l2*cos(t2)], [l1 l2*sin(t2)+l1], '-s','linewidth',1.5);
        hold on;
        text(l2*cos(t2), (l2*sin(t2)+l1), 'A2,A3');
        line([-(3*l5) (3*l5)],[l7 l7],'LineStyle','--');
        plot([0 power(-1,floor(sin(t3)))*l4*cos(t3)], [0 abs(l4*sin(t3))], '-o','linewidth',1.5);
        if(floor(sin(t3))~=0)
            t3=t3+pi;
        end
        hold on;
        text(power(-1,floor(sin(t3)))*l4*cos(t3)+l1/7, abs(l4*sin(t3)), 'B');
        t5 = asin(-(l7-l4*sin(t3))/l5);
        l6=l4*cos(t3)+l5*cos(t5);
        plot([l4*cos(t3) l6], [l4*sin(t3) l7], '-s','linewidth',1.5);
        hold on;
        text(l6+l1/7, l7, 'C');
        line([l6 l6], [l4 l2],'LineStyle','--');
%-----------------------------------------------------------------------------------------------       
        %Velocity Analysis of Crank Slotted Lever Mechanism
        subplot (2,2,2);
        w2=w0 + aa*t;
        w3= (l2*w2*cos(t2)*cos(t3) + l2*w2*sin(t2)*sin(t3))/l3;
        va2=w2*l2;
        va3=w3*l3;
        vb=abs(w3)*l4;
        w5=(-l4*w3*cos(t3))/(l5*cos(t5));
        vc=(-l4*w3*sin(t3) - l5*w5*sin(t5));
        hold off;
        quiver(0,0,1.11*va2*cos(pi/2 + t2),1.11*va2*sin(pi/2 + t2),'LineWidth',2);
        axis equal;
        axis([(-3*va2) (3*va2) (-2*va2) (2*va2)]);
        hold on;
        quiver(0,0, 10*va3*cos(pi/2 + t3)/9,10*va3*sin(pi/2 + t3)/9,'LineWidth',2);
        vsx=va3*cos(pi/2 + t3) - va2*cos(pi/2 + t2);
        vsy=va2*sin(pi/2 + t2) - va3*sin(pi/2 + t3);
        quiver(va2*cos(pi/2 + t2),va2*sin(pi/2 + t2), 1.11*vsx, -1.11*vsy,'LineWidth',2);
        text(0,0,'o');
        text(va2*cos(pi/2 + t2),va2*sin(pi/2 + t2),'a2');
        text(va3*cos(pi/2 + t3),va3*sin(pi/2 + t3),'a3');
        quiver(0,0, 1.11*vb*cos(pi/2 + t3)*va3/abs(va3), 1.11*vb*sin(pi/2 + t3)*va3/abs(va3),'LineWidth',2);
        text(vb*cos(pi/2 + t3)*va3/abs(va3),vb*sin(pi/2 + t3)*va3/abs(va3),'b');
        if cos(t2)>0
            vc=-1*vc;
        end
         quiver(0,0,1.11*vc,0,'LineWidth',2);
         text(vc,0,'c');
         vcb=sqrt( (vc)^2 + (vb)^2 - 2*vc*vb*cos(2*pi - t5) );%velocity Vcb
         plot([vb*cos(pi/2 + t3)*va3/abs(va3) vc],[vb*sin(pi/2 + t3)*va3/abs(va3) 0],'linewidth',1.5);
%-----------------------------------------------------------------------------------------------------------
        %Acceleration Analysis of Crank Slotted Lever Mechanism
        subplot (2,2,3);
        Aa2n = l2*(w2)^2; %normal acceleration of the input crank
        Aa2t = l2*aa; %tangential acceleration of the input crank
        Aa3n=l3*(w3)^2; %normal acceleration of slotted lever at point A3, on the link
        l3d= (w3*l2*cos(t2)*sin(t3) - l2*w2*cos(t3)*sin(t2))/(cos(t3)*cos(t3)); % value of l3 dot(derivative of l3)
        aa3 = (l2/l3)*( -(w2^2)*cos(t2-t3) + aa*sin(t2-t3)) - 2*l3d*w3/l3; %angluar acceleration of slotted lever
        l3dd = (1/cos(t3))*( -l2*aa*sin(t2) - l2*((w2)^2)*cos(t2) + 2*l3d*w3*sin(t3)...
            + l3*((w3)^2)*cos(t3) + l3*aa3*sin(t3)); % value of l3 double dot(slider acceleration)
        aa5 = (l5*w5*w5*sin(t5) + l4*(w3*w3*sin(t3) - aa3*cos(t3)))/(l5*cos(t5));
        Aa3t = l3*aa3;
        Abn = l4*(w3^2);
        Abt = l4*abs(aa3);
        Acbn = l5*(w5^2);
        Acbt = l5*abs(aa5);
        Ac = -l3*aa3*sin(t3) - l3*w3*w3*cos(t3) - l5*aa5*sin(t5) - l5*w5*w5*cos(t5);
        hold off;
        quiver(0,0,1.11*Aa2n*cos(-pi+t2),1.11*Aa2n*sin(-pi+t2),'LineWidth',2);
        text(0,0,'o');
        axis equal;
        axis([(-22*va2) (22*va2) (-7*va2) (7*va2)]);
        hold on;
        text(Aa2n*cos(-pi+t2), Aa2n*sin(-pi+t2), 'a2');
        quiver(Aa2n*cos(-pi+t2), Aa2n*sin(-pi+t2), 1.11*Aa2t*cos(t2+pi/2), 1.11*Aa2t*sin(t2+pi/2),'LineWidth',2);
        if t3>pi/2 && t3<3*pi/4
            Aa3n= -1*Aa3n;
        end
        as=sqrt( ((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2)) - l3dd*cos(t3))^2 + ((-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2)) - l3dd*sin(t3))^2 );
        acr= sqrt( ((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2) - l3dd*cos(t3)) - (Aa2n*cos(-pi+t2) + Aa2t*cos(t2+pi/2) ))^2 + ((-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2) -l3dd*sin(t3)) - (Aa2n*sin(-pi+t2)+ Aa2t*sin(t2+pi/2)))^2 );
        quiver(0,0,-1.11*Aa3n*cos(-pi+t3), -abs(1.11*Aa3n*sin(-pi+t3)),'LineWidth',2);%normal component of Aa3
        text(-1.11*Aa3n*cos(-pi+t3), -abs(1.11*Aa3n*sin(-pi+t3)), 'a31');
        quiver(-Aa3n*cos(-pi+t3),  -abs(Aa3n*sin(-pi+t3)), 1.11*Aa3t*cos(t3-3*pi/2), 1.11*Aa3t*sin(t3-3*pi/2),'LineWidth',2);%tangential component of Aa3
        text(-Aa3n*cos(-pi+t3)+1.11*Aa3t*cos(t3-3*pi/2), -abs(Aa3n*sin(-pi+t3))+1.11*Aa3t*sin(t3-3*pi/2), 'a32');
        quiver((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2)), (-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2)), -1.11*l3dd*cos(t3), -1.11*l3dd*sin(t3),'LineWidth',2);%slider acceleration component
        plot([(-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2) - l3dd*cos(t3)) (Aa2n*cos(-pi+t2) + Aa2t*cos(t2+pi/2) )], [(-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2) -l3dd*sin(t3)) (Aa2n*sin(-pi+t2)+ Aa2t*sin(t2+pi/2))], '-','linewidth',1.5);%coriolis acceleratoin component
        quiver(0,0, 1.11*Abn*cos(pi+t3), 1.11*Abn*sin(pi+t3),'LineWidth',2); %normal component of acceleration at point B
        text(Abn*cos(pi+t3), Abn*sin(pi+t3), 'Abn');
        quiver(Abn*cos(pi+t3), Abn*sin(pi+t3), 1.11*Abt*cos(t3-3*pi/2), 1.11*Abt*sin(t3-3*pi/2), 'LineWidth',2);
        text(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2), 'Abt');
        quiver(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2), 1.11*Acbn*cos(t5-pi), 1.11*Acbn*sin(t5-pi),'LineWidth',2);
        if cos(t2) < 0
            t5=t5-pi;
        end
        quiver(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2)+Acbn*cos(t5-pi), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2)+Acbn*sin(t5-pi), 1.11*Acbt*cos(t5-pi/2), 1.11*Acbt*sin(t5-pi/2),'LineWidth',2);
        plot([0 Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2)+Acbn*cos(t5-pi)+Acbt*cos(t5-pi/2)],[0 Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2)+Acbn*sin(t5-pi)+Acbt*sin(t5-pi/2)], '-','linewidth',1.5);
        t=t+0.05; 
 %-----------------------------------------------------------------------------------------%
        %Dynamic force analysis of Crank slotted lever mechanism%
            %normal simulation of assembly%
            subplot (2,2,4);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold off;
            plot([0 0],[0 l1], '--o','linewidth',1.5);
            axis equal;
            axis([-(3*l5) (3*l5) -l1/2 (l1+l4)]);
            t2=(w0*t + 0.5*aa*(t*t));
            t3=pi+atan( (l2*sin(t2)+l1)/(l2*cos(t2)) );
            l3=l2 * cos(t2) / cos(t3);
            plot([0 l2*cos(t2)], [l1 l2*sin(t2)+l1], '-s','linewidth',1.5);
            hold on;
            line([-(3*l5) (3*l5)],[l7 l7],'LineStyle','--');
            plot([0 power(-1,floor(sin(t3)))*l4*cos(t3)], [0 abs(l4*sin(t3))], '-o','linewidth',1.5);
            if(floor(sin(t3))~=0)
                t3=t3+pi;
            end
            hold on;
            t5 = asin(-(l7-l4*sin(t3))/l5);
            l6=l4*cos(t3)+l5*cos(t5);
            plot([l4*cos(t3) l6], [l4*sin(t3) l7], '-s','linewidth',1.5);
            hold on;
            line([l6 l6], [l4 l2],'LineStyle','--');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            quiver(l6+l2,l7, -l2, 0, 'linewidth', 3.5);
            Ga2n=Aa2n/2;
            Ga2t=Aa2t/2;
            Ga2=sqrt(Ga2n^2+Ga2t^2);
            Ga3n=Aa3n/2;
            Ga3t=Aa3t/2;
            Ga3=sqrt(Ga3n^2+Ga3t^2);
            Gbn=Abn/2;
            Gbt=Abt/2;
            Gb=sqrt(Gbn^2+Gbt^2);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %from the FBD of slider%
%             F56n*cos(2*pi-t5) - F - F56t*cos(t5-pii/2)= m*Ac;
%             N-F56n*sin(2*pi-t5)-F56t*sin(t5-pi/2)=0;
%             %from the FBD of link 5%
%             F65n=-F56n;
%             F65t=-F56t;
              I5=(rho*l5^3)/12;
%             F45n - F65n=(rho*l5)*Gbn;
%             -F45t-F65t=Gbt*(rho*l5);
%             -F45t*l5/2 + F65t*l5/2 = I5*aa5;
%             %from the FBD of link 4%
%             F54n=-F45n;
%             F54t=-F45n;
%             F54t+F14n=(rho*l4)*Ga3n;
%             F54n-F14t-F34=(rho*l4)*Ga3t;
%             (F54n+F14-F34)*l4/2 + F34*l3=((rho*l4^3)/12)*aa3;
%             %fbd of link 3%
%             F43=-F34;
%             F23n*cos(t3-t2) + F23t*cos(pi/2-t3+t2) = m*as;
%             F43+F23t*sin(pi/2-t3+t2)-F23n*sin(t3-t2)=m*acr;
%             %fbd of the crank%
%             F32n=-F23n;
%             F32t=-F23t;
%             F12n-F32n=(rho*l2)*Ga2n;
%             -F12t-F32t=(rho*l2)*Ga2t;
%             T+(F12t-F32t)*l2/2=((rho*l2^3)/12)*aa2;
%initialization of variables%
F56n= 0;F56t=0; N=0; F65n=0; F65t=0; F45n=0; F45t=0; F54n=0; 
F54t=0; F14n=0; F14t=0; F34=0; F43=0; F23n=0; F23t=0; F32n=0; F32t=0; F12n=0; F12t=0; T=0;
            % 20 x 20 matrix comprising above equations to be solved%
            A=[cos(2*pi-t5) -cos(t5-pi/2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               -sin(2*pi-t5) -sin(t5-pi/2) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               1              0             0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               0              1             0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               0              0             0 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 -l5/2 0 l5/2 0 0 0 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 0 1 0 0 -1 0 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 0 l4/2 0 0 l4/2 (l3-l4/2) 0 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 cos(t3-t2) cos(pi/2-t3+t2) 0 0 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 1 -sin(t3-t2) sin(pi/2-t3+t2) 0 0 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 0;...
               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (-l2/2) 0 l2/2 1];
           X=[F56n; F56t; N; F65n; F65t; F45n; F45t; F54n; F54t; F14n; F14t; F34; F43; F23n; F23t; F32n; F32t; F12n; F12t; T];
           B=[(F+m*Ac); 0; 0; 0; (rho*l5*aa5); (Gbt*rho*l5); (I5*aa5); 0;0; (rho*l4*Ga3n); (rho*l4*Ga3t); (rho*l4^3*aa3)/12; 0; m*as; m*acr; 0; 0; (rho*l2*Ga2n); (rho*l2*Ga2t); (rho*l2^3*aa)/12];
           X=inv(A)*B;
           fprintf('Torque required = %f\n Shaking Force at O2, F21n=%f and F21t=%f\nNOTE:F21n and F21t are in the normal and tangential direction to the input link\n',X(20,1),-X(18,1),-X(19,1));
    end    
    close %closes the plot window
end
elseif choice==4
    %Taking input from the user and storing in variables
    l1=str2double(ss(1));
    l2=str2double(ss(2));
    l4=str2double(ss(3));
    l5=str2double(ss(4));
    l7=str2double(ss(5));
    t=str2double(ss(6));
    w0=str2double(ss(7));
    aa=str2double(ss(8));
    F=str2double(ss(9));
    %checking validity of link lengths
    if (l4 < 2*l2) || (((l4 - l5) > l7) || (l5 > l4)) || (l1 + l2 > l4) || (l7 > 0.707*l4 + l5) || (l2>l1) || l2<=0 || l1<=0 || l4<=0 || l5<=0 || l7<=0
    h = errordlg('The links are not valid to form a mechanism.','Error');
    return
    else
    h=figure;
        %Position Analysis and Configuration Diagram of Crank Slotted Lever Mechanism%
        subplot(2,2,1);
        text(l1/10,0, 'O2');
        plot([0 0],[0 l1], '--o','linewidth',1.5);
        hold on;
        text(-l1/3,l1, 'O4');
        axis equal;
        axis([-(3*l5) (3*l5) -l1/2 (l1+l4)]);
        t2=t*pi/180;
        t3=pi+atan( (l2*sin(t2)+l1)/(l2*cos(t2)) );
        l3=l2 * cos(t2) / cos(t3);
        plot([0 l2*cos(t2)], [l1 l2*sin(t2)+l1], '-s','linewidth',1.5);
        hold on;
        text(l2*cos(t2), (l2*sin(t2)+l1), 'A2,A3');
        line([-(3*l5) (3*l5)],[l7 l7],'LineStyle','--');
        plot([0 power(-1,floor(sin(t3)))*l4*cos(t3)], [0 abs(l4*sin(t3))], '-o','linewidth',1.5);
        if(floor(sin(t3))~=0)
            t3=t3+pi;
        end
        hold on;
        text(power(-1,floor(sin(t3)))*l4*cos(t3)+l1/7, abs(l4*sin(t3)), 'B');
        t5 = asin(-(l7-l4*sin(t3))/l5);
        l6=l4*cos(t3)+l5*cos(t5);
        plot([l4*cos(t3) l6], [l4*sin(t3) l7], '-s','linewidth',1.5);
        hold on;
        text(l6+l1/7, l7, 'C');
        line([l6 l6], [l4 l2],'LineStyle','--');
%-----------------------------------------------------------------------------------------------       
        %Velocity Analysis of Crank Slotted Lever Mechanism%
        subplot (2,2,2);
        w2=w0;
        w3= (l2*w2*cos(t2)*cos(t3) + l2*w2*sin(t2)*sin(t3))/l3;
        va2=w2*l2;
        va3=w3*l3;
        vb=abs(w3)*l4;
        w5=(-l4*w3*cos(t3))/(l5*cos(t5));
        vc=(-l4*w3*sin(t3) - l5*w5*sin(t5));
        hold off;
        quiver(0,0,1.11*va2*cos(pi/2 + t2),1.11*va2*sin(pi/2 + t2),'LineWidth',2);
        axis equal;
        axis([(-3*va2) (3*va2) (-2*va2) (2*va2)]);
        hold on;
        quiver(0,0, 10*va3*cos(pi/2 + t3)/9,10*va3*sin(pi/2 + t3)/9,'LineWidth',2);
        vsx=va3*cos(pi/2 + t3) - va2*cos(pi/2 + t2);
        vsy=va2*sin(pi/2 + t2) - va3*sin(pi/2 + t3);
        quiver(va2*cos(pi/2 + t2),va2*sin(pi/2 + t2), 1.11*vsx, -1.11*vsy,'LineWidth',2);
        text(0,0,'o');
        text(va2*cos(pi/2 + t2),va2*sin(pi/2 + t2),'a2');
        text(va3*cos(pi/2 + t3),va3*sin(pi/2 + t3),'a3');
        quiver(0,0, 1.11*vb*cos(pi/2 + t3)*(va3/abs(va3)), 1.11*vb*sin(pi/2 + t3)*(va3/abs(va3)),'LineWidth',2);
        text(vb*cos(pi/2 + t3)*va3/abs(va3),vb*sin(pi/2 + t3)*va3/abs(va3),'b');
        if cos(t2)>0
            vc=-1*vc;
         end
         quiver(0,0,1.11*vc,0,'LineWidth',2);
         text(vc,0,'c');
         vcb=sqrt( (vc)^2 + (vb)^2 - 2*vc*vb*cos(2*pi - t5) ); %velocity Vcb
         plot([vb*cos(pi/2 + t3)*va3/abs(va3) vc],[vb*sin(pi/2 + t3)*va3/abs(va3) 0],'linewidth',1.5);
%-----------------------------------------------------------------------------------------------------------
        %Acceleration Analysis of Crank Slotted Lever Mechanism
        subplot (2,2,3);
        Aa2n = l2*(w2)^2; %normal acceleration of the input crank
        Aa2t = l2*aa; %tangential acceleration of the input crank
        Aa3n=l3*(w3)^2; %normal acceleration of slotted lever at point A3, on the link
        l3d= (w3*l2*cos(t2)*sin(t3) - l2*w2*cos(t3)*sin(t2))/(cos(t3)*cos(t3)); % value of l3 dot(derivative of l3)
        aa3 = (l2/l3)*( -(w2^2)*cos(t2-t3) + aa*sin(t2-t3)) - 2*l3d*w3/l3; %angluar acceleration of slotted lever
        l3dd = (1/cos(t3))*( -l2*aa*sin(t2) - l2*((w2)^2)*cos(t2) + 2*l3d*w3*sin(t3) + l3*((w3)^2)*cos(t3) + l3*aa3*sin(t3)); % value of l3 double dot(slider acceleration)
        aa5 = (l5*w5*w5*sin(t5) + l4*(w3*w3*sin(t3) - aa3*cos(t3)))/(l5*cos(t5));
        Aa3t = l3*aa3;
        Abn = l4*(w3^2);
        Abt = l4*abs(aa3);
        Acbn = l5*(w5^2);
        Acbt = l5*abs(aa5);
        Ac = -l3*aa3*sin(t3) - l3*w3*w3*cos(t3) - l5*aa5*sin(t5) - l5*w5*w5*cos(t5);
        hold off;
        quiver(0,0,1.11*Aa2n*cos(-pi+t2),1.11*Aa2n*sin(-pi+t2),'LineWidth',2);
        axis equal;
        axis([(-13*va2) (13*va2) (-8*va2) (8*va2)]);
        hold on;
        text(Aa2n*cos(-pi+t2)/2, Aa2n*sin(-pi+t2)/2, 'Aa2n');
        quiver(Aa2n*cos(-pi+t2), Aa2n*sin(-pi+t2), 1.11*Aa2t*cos(t2+pi/2), 1.11*Aa2t*sin(t2+pi/2),'LineWidth',2);
        if t3>pi/2 && t3<3*pi/4
            Aa3n= -1*Aa3n;
        end
        quiver(0,0,-1.11*Aa3n*cos(-pi+t3), -abs(1.11*Aa3n*sin(-pi+t3)),'LineWidth',2);
        text(0,0,'o');
        text(-1.11*Aa3n*cos(-pi+t3), -abs(1.11*Aa3n*sin(-pi+t3)), 'a31');
        quiver(-Aa3n*cos(-pi+t3),  -abs(Aa3n*sin(-pi+t3)), 1.11*Aa3t*cos(t3-3*pi/2), 1.11*Aa3t*sin(t3-3*pi/2),'LineWidth',2);
        text(-Aa3n*cos(-pi+t3)+1.11*Aa3t*cos(t3-3*pi/2), -abs(Aa3n*sin(-pi+t3))+1.11*Aa3t*sin(t3-3*pi/2), 'a32');
        quiver((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2)), (-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2)), -1.11*l3dd*cos(t3), -1.11*l3dd*sin(t3),'LineWidth',2);
         quiver(0,0, 1.11*Abn*cos(pi+t3), 1.11*Abn*sin(pi+t3),'LineWidth',2); %normal component of acceleration at point B
        text(Abn*cos(pi+t3), Abn*sin(pi+t3), 'Abn');
        quiver(Abn*cos(pi+t3), Abn*sin(pi+t3), 1.11*Abt*cos(t3-3*pi/2), 1.11*Abt*sin(t3-3*pi/2), 'LineWidth',2);
        text(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2), 'Abt');
        quiver(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2), 1.11*Acbn*cos(t5-pi), 1.11*Acbn*sin(t5-pi),'LineWidth',2);
        if cos(t2) < 0
            t5=t5-pi;
        end
        quiver(Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2)+Acbn*cos(t5-pi), Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2)+Acbn*sin(t5-pi), 1.11*Acbt*cos(t5-pi/2), 1.11*Acbt*sin(t5-pi/2),'LineWidth',2);
        plot([0 Abn*cos(pi+t3) + Abt*cos(t3-3*pi/2)+Acbn*cos(t5-pi)+Acbt*cos(t5-pi/2)],[0 Abn*sin(pi+t3)+Abt*sin(t3-3*pi/2)+Acbn*sin(t5-pi)+Acbt*sin(t5-pi/2)], '-','linewidth',1.5);
        as=sqrt( ((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2)) - l3dd*cos(t3))^2 + ((-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2)) - l3dd*sin(t3))^2 );
        plot([(-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2) - l3dd*cos(t3)) (Aa2n*cos(-pi+t2) + Aa2t*cos(t2+pi/2) )], [(-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2) -l3dd*sin(t3)) (Aa2n*sin(-pi+t2)+ Aa2t*sin(t2+pi/2))], '-','linewidth',1.5);  
        acr= sqrt( ((-Aa3n*cos(-pi+t3) + Aa3t*cos(t3-3*pi/2) - l3dd*cos(t3)) - (Aa2n*cos(-pi+t2) + Aa2t*cos(t2+pi/2) ))^2 + ((-abs(Aa3n*sin(-pi+t3)) + Aa3t*sin(t3-3*pi/2) -l3dd*sin(t3)) - (Aa2n*sin(-pi+t2)+ Aa2t*sin(t2+pi/2)))^2 );
        %Creating the tabular column
        f=msgbox('Maximize the plot window for better visibility; Please ZOOM and PAN into the plot to view vectors clearly', 'Info', 'help');
        cname={'Magnitude'};
        rname={'Va2'; 'Va3'; 'Va3a2'; 'Vb'; 'Vc'; 'Vcb';'Aa3n'; 'Aa3t';'Aa2n';'Aa2t';'Aa2a3CR';'Aa2a3S';'Abn';'Abt';'Acbn';'Acbt';'Ac'};
        d=[abs(va2);abs(va3);sqrt(vsx^2 + vsy^2); abs(vb); abs(vc); abs(vcb); abs(Aa3n);abs(Aa3t);abs(Aa2n);abs(Aa2t); abs(acr); abs(as);Abn;Abt;Acbn;Acbt;abs(Ac)];
        z=uitable(h,'Data',d, 'ColumnName',cname,'RowName',rname);
        z.Position(3) = z.Extent(3);
        z.Position(4) = z.Extent(4);
        fgcolor = z.ForegroundColor;
        z.ForegroundColor = [1 0.25 0.56];
%---------------------------------------------------------------------------%
        %Static Force Analysis of a crank slotted mechanism%
        subplot 224;
        plot([0 0],[0 l1], '--o','linewidth',1.5);
        hold on;
        axis equal;
        axis([-(3*l5) (3*l5) -l1/2 (l1+l4)]);
        t2=t*pi/180;
        t3=pi+atan( (l2*sin(t2)+l1)/(l2*cos(t2)) );
        l3=l2 * cos(t2) / cos(t3);
        plot([0 l2*cos(t2)], [l1 l2*sin(t2)+l1], '-s','linewidth',1.5);
        hold on;
        line([-(3*l5) (3*l5)],[l7 l7],'LineStyle','--');
        plot([0 power(-1,floor(sin(t3)))*l4*cos(t3)], [0 abs(l4*sin(t3))], '-o','linewidth',1.5);
        if(floor(sin(t3))~=0)
            t3=t3+pi;
        end
        hold on;
        t5 = asin(-(l7-l4*sin(t3))/l5);
        l6=l4*cos(t3)+l5*cos(t5);
        plot([l4*cos(t3) l6], [l4*sin(t3) l7], '-s','linewidth',1.5);
        hold on;
        line([l6 l6], [l4 l2],'LineStyle','--');
        quiver(l6+l2,l7, -l2, 0, 'linewidth', 3.5);
        text(l6+l2/2, l7+l2/3, 'F');
        %actual force analysis starts
        F56=F/cos(t5);
        F54=-F56;
        M=[-sin(t3) 1 0; cos(t3) 0 1; (l3-l4) -l4*sin(t3) l4*cos(t3)];
        N=[ F54*cos(t5); -F54*sin(t5); 0];
        O=inv(M)*N;
        F34=O(1,1);
        F14x=O(2,1);
        F14y=O(3,1);
        F32=-F34;
        T=F32*sin(pi/2 - t2 + t3); %required torque%
        if T<0
            fprintf('The value of TORQUE is %f CLOCKWISE\n',abs(T));
        else
            fprintf('The value of TORQUE is %f ANTICLOCKWISE\n',abs(T));
        end
        f=msgbox('Please view the command window for the REQUIRED TORQUE VALUE and DIRECTION');     
    end
else
    return
end
%End of the code%