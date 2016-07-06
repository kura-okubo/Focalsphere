% Focalsphere(azimuth,dip,rake) plots focal sphere and projections.
% Syntax:   Focalsphere(azimuth,dip,rake) or Focalsphere()
%
% Inputs:   azimuth, dip and rake for the focal sphere.
%           Each unit should be degree. If there is no input, it plots sample case.
% 
% Notation:
%           n   normal unit vector of the fault plane
%           I   Slip direction
%           z   Pole of the focal shere
%           red and blue plane show the fault and auxiliary plane
%           respectively.
%
% You can play three parameters (azimuth, dip and rake) with sliders. 
%
% Reference:
% Udias, A., Madariaga, R., Buforn, E.(2014), Source Mechanisms of Earthquakes: Theory and Practice, Cambridge University Press. 
%
% Copyright (c) 2016 Kurama OKUBO.
% $Date: 2016/04$

function  Focalsphere(azimuth,dip,rake)

if nargin == 3
    
    %initial parameters
    if azimuth < 0 || azimuth > 360 ||...
            dip <-90 || dip> 90 ||...
            rake < 0 || rake > 360
        error('Value is outside of range. please set within 0< azimuth&rake <360 and -90< dip <90.');                
    end
    %%%%%%%%%%%%%%%
    phi = azimuth; %azimuth
    del = dip; %dip
    lamb = rake; %rake
    %%%%%%%%%%%%%%%
else
     %initial parameters
    %%%%%%%%%%%%%%%
    phi = 30; %azimuth
    del = 50; %dip
    lamb = 20; %rake
    %%%%%%%%%%%%%%%
end

%%%color vector%%
re = [1.000 ,0.2627 ,0.000];
bl = [0.000 ,0.3176 ,1.000];
ye = [1.000 ,0.8627 ,0.000];

%maniplate angles
phi = -deg2rad(phi);
del = deg2rad(del);
lamb = deg2rad(lamb);

%*****************************************************%
%Unit vector and fault and axiliary plane
%*****************************************************%

%make unit of I and normal N
I = [cos(lamb)*cos(phi)+cos(del)*sin(lamb)*sin(phi); cos(lamb)*sin(phi)-cos(del)*sin(lamb)*cos(phi); -sin(lamb)*sin(del)];
N = [sin(del)*sin(phi); -sin(del)*cos(phi); cos(del)];
zp = cross(I,N);

%make fault plane
lamb0=deg2rad((0:1:360));
I0 = [cos(lamb0).*cos(phi)+cos(del).*sin(lamb0)*sin(phi); cos(lamb0).*sin(phi)-cos(del).*sin(lamb0)*cos(phi); -sin(lamb0).*sin(del)];
p = [0, 0, 0];

%make auxiliary plane
for i = 1:size(I0,2)
    x0 = I0(1,i);
    y0 = I0(2,i);
    z0 = I0(3,i);
    [x00 y00 z00] = rotatexyz(x0,y0,z0,phi,'z');
    [x1 y1 z1] = rotatexyz(x00,y00,z00,del,'x');
    [x2 y2 z2] = rotatexyz(x1,y1,z1,pi/2-lamb,'z');
    [x3 y3 z3] = rotatexyz(x2,y2,z2,pi/2,'x');
    [x4 y4 z4] = rotatexyz(x3,y3,z3,-(pi/2-lamb),'z');
    [x5 y5 z5] = rotatexyz(x4,y4,z4,-del,'x');
    [Aux(1,i) Aux(2,i) Aux(3,i)] = rotatexyz(x5,y5,z5,-phi,'z');
end

%*****************************************************%
%Horizontal Stereographic projection of the lower hemisphere of the focal sphere
%*****************************************************%

%extract lower hemisphere
I1 = I0;
Aux1 = Aux;
for i = 1:length(lamb0);
    if I1(3,i) > 0
        I1(:,i) = nan;
    end
    
    if Aux1(3,i) > 0
        Aux1(:,i) = nan;
    end
end

%*****************************************************%
%Vertical Stereographic projection of the half hemisphere of the focal sphere
%*****************************************************%

%extract lower hemisphere
I2 = I0;
Aux2 = Aux;
for i = 1:length(lamb0);
    if I2(1,i) < 0
        I2(:,i) = nan;
    end
    
    if Aux2(1,i) < 0
        Aux2(:,i) = nan;
    end
end


%*****************************************************%
%Dilatation and Compression coloring
%*****************************************************%

%make quarter sphere
[X,Y,Z] = sphere;
%You can abailable hemisphere with using under 11
hemi = 11;
icount=0;
for i = 1:hemi
    for j = 1:hemi
        icount = icount+1;
        X0(icount) = X(i,j);
        Y0(icount) = Y(i,j);
        Z0(icount) = Z(i,j);
    end
end

%rotate

for i = 1:icount
    %ajust coordinate
    [xx1 yy1 zz1] = rotatexyz(X0(i),Y0(i),Z0(i),pi/2,'x');
    [xx2 yy2 zz2] = rotatexyz(X0(i),Y0(i),Z0(i),-pi/2,'x');
    %ajust coordinate
    [xx1 yy1 zz1] = rotatexyz(xx1,yy1,zz1,-(pi/2-lamb),'z');
    [xx1 yy1 zz1] = rotatexyz(xx1,yy1,zz1,-del,'x');
    [X1(i) Y1(i) Z1(i)] = rotatexyz(xx1,yy1,zz1,-phi,'z');
    [xx2 yy2 zz2] = rotatexyz(xx2,yy2,zz2,-(pi/2-lamb),'z');
    [xx2 yy2 zz2] = rotatexyz(xx2,yy2,zz2,-del,'x');
    [X2(i),Y2(i),Z2(i)] = rotatexyz(xx2,yy2,zz2,-phi,'z');
end

%revert to mesh
icount = 0;
for i = 1:hemi
    for j = 1:hemi
        icount = icount+1;
        XR(i,j)=X1(icount);
        YR(i,j)=Y1(icount);
        ZR(i,j)=Z1(icount);
        XL(i,j)=X2(icount);
        YL(i,j)=Y2(icount);
        ZL(i,j)=Z2(icount);
    end
end

%*****************************************************%
%plot figure
%*****************************************************%

f = figure(1);
clf;


hold on;
f.Position(3) = 650; 
view([-60 30]);
xlim([-1 1.5]);
ylim([-1 1]);
zlim([-1.5 1]);
ax = gca;
ax.XTick = [-1 0 1];
ax.XTickLabel = {};
ax.YTick = [-1 0 1];
ax.YTickLabel = {};
ax.ZTick = [-1 0 1];
ax.ZTickLabel = {};

axis equal;

%Unit Vectors
plot3([p(1) I(1)],[p(2) I(2)],[p(3) I(3)],'Color',ye,'LineWidth',3);
plot3([p(1) N(1)],[p(2) N(2)],[p(3) N(3)],'Color',re,'LineWidth',3);
plot3([p(1) zp(1)],[p(2) zp(2)],[p(3) zp(3)],'Color',bl,'LineWidth',3);

%Fault and auxility planes
fill3(I0(1,:),I0(2,:),I0(3,:),re, 'FaceAlpha',.5);
fill3(Aux(1,:),Aux(2,:),Aux(3,:),bl,'FaceAlpha',.5);

%sphere grid
[X,Y,Z] = sphere;
surf(X,Y,Z,'FaceColor','none','EdgeColor','k','EdgeAlpha',0.1);
grid on;

%Horizontal projection
xc = cos(lamb0);
yc = sin(lamb0);
zc = -1.5*ones(length(lamb0),1);
plot3(xc,yc,zc,'k-');
plot3(I1(1,:),I1(2,:),zc,'Color',re);
plot3(Aux1(1,:),Aux1(2,:),zc,'Color',bl);

%Vertical projection
xc = 1.5*ones(length(lamb0),1) ;
yc = cos(lamb0);
zc = sin(lamb0);
plot3(xc,yc,zc,'k-');
plot3(xc,I2(2,:),I2(3,:),'Color',re);
plot3(xc,Aux2(2,:),Aux2(3,:),'Color',bl);

%refference axis
plot3([-1 1], [0 0],[0 0],'k--');
plot3([0 0], [-1 1],[0 0],'k--');
plot3([0 0], [0 0],[-1 1],'k--');


%black color for focal sphere
surf(XR,YR,ZR,'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none');
surf(XL,YL,ZL,'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none');


%Texts
sazim = sprintf('azimuth=%4.1f', -rad2deg(phi));
sdip = sprintf('dip=%4.1f', rad2deg(del));
srake = sprintf('rake=%4.1f', rad2deg(lamb));
title([sazim '$^\circ$ ' sdip '$^\circ$ ' srake '$^\circ$'],'Interpreter','LaTex', 'FontSize', 20);
xlabel('NS');
ylabel('EW');
zlabel('Z');
text(I(1),I(2),I(3),'$\mathbf{I}$','Interpreter','LaTex','FontSize',24)
text(N(1),N(2),N(3),'$\mathbf{n}$','Interpreter','LaTex','FontWeight','bold','FontSize',24)
text(zp(1),zp(2),zp(3),'$\mathbf{z}$','Interpreter','LaTex','FontWeight','bold','FontSize',24)
text(1,0,0,'N','FontSize',24)
text(-1,0,0,'S','FontSize',24)
text(0,1,0,'W','FontSize',24)
text(0,-1,0,'E','FontSize',24)
text(0,0,1,'Z','FontSize',24)


%*****************************************************%
%uicontrol
%*****************************************************%

f.Visible = 'off';

setappdata(f,'m_phi',-rad2deg(phi));
setappdata(f,'m_del',rad2deg(del));
setappdata(f,'m_lamb',rad2deg(lamb));

%azimuth slider
s_azimuth = uicontrol('Parent', f, 'Style', 'slider',...
        'Min',0,'Max',360,'Value',rad2deg(-phi),...
        'Tag','azimuth',...
        'Position', [20 180 140 20],...
        'Callback', {@FScontrol,'azimuth'}); 
txt = uicontrol('Style','text',...
        'FontSize',16,...
        'Position',[20 200 140 20],...
        'String','azimuth');

%dip slider
s_dip = uicontrol('Parent', f,'Style', 'slider',...
        'Min',-90,'Max',90,'Value',rad2deg(del),...
        'Tag','dip',...
        'Position', [20 140 140 20],...
        'Callback', {@FScontrol,'dip'});
txt = uicontrol('Style','text',...
        'FontSize',16,...
        'Position', [20 160 140 20],...
        'String','dip');

%rake slider
s_rake = uicontrol('Parent', f,'Style', 'slider',...
        'Min',0,'Max',360,'Value',rad2deg(lamb),...
        'Tag','rake',...
        'Position',[20 100 140 20],...
        'Callback', {@FScontrol,'rake'});

txt = uicontrol('Style','text',...
        'Position',[20 120 140 20],...
        'FontSize',16,...
        'String','rake');
    
f.Visible = 'on';

    function FScontrol(source,callbackdata,type)      
        
        switch char(type)
            case 'azimuth'
                    phi = source.Value;
                    del = getappdata(f,'m_del');
                    lamb = getappdata(f,'m_lamb');
                    setappdata(f,'m_phi',source.Value);
            case  'dip'
                    phi = getappdata(f,'m_phi');
                    del = source.Value;
                    lamb = getappdata(f,'m_lamb');
                    setappdata(f,'m_del',source.Value);

            case  'rake'
                    phi = getappdata(f,'m_phi');
                    del = getappdata(f,'m_del');
                    lamb = source.Value;
                    setappdata(f,'m_lamb',source.Value);
        end  

        
    
        % maniplate angl%return modified vale  
        
        
        phi = -deg2rad(phi);
        del = deg2rad(del);
        lamb = deg2rad(lamb);

        %make unit of I and normal N
        I = [cos(lamb)*cos(phi)+cos(del)*sin(lamb)*sin(phi); cos(lamb)*sin(phi)-cos(del)*sin(lamb)*cos(phi); -sin(lamb)*sin(del)];
        N = [sin(del)*sin(phi); -sin(del)*cos(phi); cos(del)];
        zp = cross(I,N);

        %make fault plane
        lamb0=deg2rad((0:1:360));
        I0 = [cos(lamb0).*cos(phi)+cos(del).*sin(lamb0)*sin(phi); cos(lamb0).*sin(phi)-cos(del).*sin(lamb0)*cos(phi); -sin(lamb0).*sin(del)];
        p = [0, 0, 0];

        %make auxiliary plane
        for i = 1:size(I0,2)
            x0 = I0(1,i);
            y0 = I0(2,i);
            z0 = I0(3,i);
            [x00 y00 z00] = rotatexyz(x0,y0,z0,phi,'z');
            [x1 y1 z1] = rotatexyz(x00,y00,z00,del,'x');
            [x2 y2 z2] = rotatexyz(x1,y1,z1,pi/2-lamb,'z');
            [x3 y3 z3] = rotatexyz(x2,y2,z2,pi/2,'x');
            [x4 y4 z4] = rotatexyz(x3,y3,z3,-(pi/2-lamb),'z');
            [x5 y5 z5] = rotatexyz(x4,y4,z4,-del,'x');
            [Aux(1,i) Aux(2,i) Aux(3,i)] = rotatexyz(x5,y5,z5,-phi,'z');
        end

        %extract lower hemisphere
        I1 = I0;
        Aux1 = Aux;
        for i = 1:length(lamb0);
            if I1(3,i) > 0
                I1(:,i) = nan;
            end

            if Aux1(3,i) > 0
                Aux1(:,i) = nan;
            end
        end

        %extract lower hemisphere
        I2 = I0;
        Aux2 = Aux;
        for i = 1:length(lamb0);
            if I2(1,i) < 0
                I2(:,i) = nan;
            end

            if Aux2(1,i) < 0
                Aux2(:,i) = nan;
            end
        end

        %make quarter sphere
        [X,Y,Z] = sphere;
        %You can abailable hemisphere with using under 11
        hemi = 11;
        icount=0;
        for i = 1:hemi
            for j = 1:hemi
                icount = icount+1;
                X0(icount) = X(i,j);
                Y0(icount) = Y(i,j);
                Z0(icount) = Z(i,j);
            end
        end

        %rotate

        for i = 1:icount
            %ajust coordinate
            [xx1 yy1 zz1] = rotatexyz(X0(i),Y0(i),Z0(i),pi/2,'x');
            [xx2 yy2 zz2] = rotatexyz(X0(i),Y0(i),Z0(i),-pi/2,'x');
            %ajust coordinate
            [xx1 yy1 zz1] = rotatexyz(xx1,yy1,zz1,-(pi/2-lamb),'z');
            [xx1 yy1 zz1] = rotatexyz(xx1,yy1,zz1,-del,'x');
            [X1(i) Y1(i) Z1(i)] = rotatexyz(xx1,yy1,zz1,-phi,'z');
            [xx2 yy2 zz2] = rotatexyz(xx2,yy2,zz2,-(pi/2-lamb),'z');
            [xx2 yy2 zz2] = rotatexyz(xx2,yy2,zz2,-del,'x');
            [X2(i),Y2(i),Z2(i)] = rotatexyz(xx2,yy2,zz2,-phi,'z');
        end

        %revert to mesh
        icount = 0;
        for i = 1:hemi
            for j = 1:hemi
                icount = icount+1;
                XR(i,j)=X1(icount);
                YR(i,j)=Y1(icount);
                ZR(i,j)=Z1(icount);
                XL(i,j)=X2(icount);
                YL(i,j)=Y2(icount);
                ZL(i,j)=Z2(icount);
            end
        end

        %*****************************************************%
        %plot figure
        %*****************************************************%
        cla (gca);
        %Unit Vectors
        plot3([p(1) I(1)],[p(2) I(2)],[p(3) I(3)],'Color',ye,'LineWidth',3);
        plot3([p(1) N(1)],[p(2) N(2)],[p(3) N(3)],'Color',re,'LineWidth',3);
        plot3([p(1) zp(1)],[p(2) zp(2)],[p(3) zp(3)],'Color',bl,'LineWidth',3);

        %Fault and auxility planes
        fill3(I0(1,:),I0(2,:),I0(3,:),re, 'FaceAlpha',.5);
        fill3(Aux(1,:),Aux(2,:),Aux(3,:),bl,'FaceAlpha',.5);

        %sphere grid
        [X,Y,Z] = sphere;
        surf(X,Y,Z,'FaceColor','none','EdgeColor','k','EdgeAlpha',0.1);
        grid on;

        %Horizontal projection
        xc = cos(lamb0);
        yc = sin(lamb0);
        zc = -1.5*ones(length(lamb0),1);
        plot3(xc,yc,zc,'k-');
        plot3(I1(1,:),I1(2,:),zc,'Color',re);
        plot3(Aux1(1,:),Aux1(2,:),zc,'Color',bl);

        %Vertical projection
        xc = 1.5*ones(length(lamb0),1) ;
        yc = cos(lamb0);
        zc = sin(lamb0);
        plot3(xc,yc,zc,'k-');
        plot3(xc,I2(2,:),I2(3,:),'Color',re);
        plot3(xc,Aux2(2,:),Aux2(3,:),'Color',bl);

        %refference axis
        plot3([-1 1], [0 0],[0 0],'k--');
        plot3([0 0], [-1 1],[0 0],'k--');
        plot3([0 0], [0 0],[-1 1],'k--');


        %black color for focal sphere
        surf(XR,YR,ZR,'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none');
        surf(XL,YL,ZL,'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none');


        %Texts
        sazim = sprintf('azimuth=%4.1f', -rad2deg(phi));
        sdip = sprintf('dip=%4.1f', rad2deg(del));
        srake = sprintf('rake=%4.1f', rad2deg(lamb));
        title([sazim '$^\circ$ ' sdip '$^\circ$ ' srake '$^\circ$'],'Interpreter','LaTex', 'FontSize', 20);
        xlabel('NS');
        ylabel('EW');
        zlabel('Z');
        text(I(1),I(2),I(3),'$\mathbf{I}$','Interpreter','LaTex','FontSize',24)
        text(N(1),N(2),N(3),'$\mathbf{n}$','Interpreter','LaTex','FontWeight','bold','FontSize',24)
        text(zp(1),zp(2),zp(3),'$\mathbf{z}$','Interpreter','LaTex','FontWeight','bold','FontSize',24)
        text(1,0,0,'N','FontSize',24)
        text(-1,0,0,'S','FontSize',24)
        text(0,1,0,'W','FontSize',24)
        text(0,-1,0,'E','FontSize',24)
        text(0,0,1,'Z','FontSize',24)

    end

    function [X,Y,Z] = rotatexyz(x,y,z,theta,mode);

        s = sin(theta);
        c = cos(theta);

        switch char(mode)    
            case 'x'
                R = [1, 0 , 0;
                    0, c, s;
                    0, -s, c];
                A = [x, y, z];
                temp = R*A';
                X = temp(1);
                Y = temp(2);
                Z = temp(3);   

            case 'y'
                R = [c, 0 , -s;
                    0, 1, 0;
                    s, 0, c];
                B = [x, y, z];
                temp = R*B';
                X = temp(1);
                Y = temp(2);
                Z = temp(3);    
            case 'z'
                R = [c, s , 0;
                    -s, c, 0;
                    0, 0, 1];
                C = [x, y, z];
                temp = R*C';
                X = temp(1);
                Y = temp(2);
                Z = temp(3);

            otherwise
                disp('error at rotatexyz; please input mode');
                return;
        end
    end
end
