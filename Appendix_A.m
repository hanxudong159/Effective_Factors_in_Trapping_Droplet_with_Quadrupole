clear all;close all; clc
N = 30;                     % Order of the Hill's Determinant
T = cell(1,4);              % Cell array to store Hill determinant matrices
EE = cell(1,4);
epsilon = -20:0.05:20;
for i = 1:length(epsilon)
    e = epsilon(i);
    a = e*ones(N,1);        % Sub diagonal elements
    c = e*ones(N,1);        % Super diagonal elements
    % Hill determinants for odd sine and cosine coefficients
    ndeto = 2*N+1;
    bo = (1:2:ndeto).^2;    % Diagonal elements
    T1 = diag(a,-1)+diag(bo)+diag(c,+1);
    T1(1,1) = 1+e;
    T{1} = T1 ;             % Odd cosine determinant matrix
    T1(1,1) = 1-e;
    T2 = T1;
    T{2} = T2 ;             % Odd sine determinant matrix
    % Hill determinants for even sine and cosine coefficients
    ndete = 2*(N+1);
    be = (2:2:ndete).^2;    % Diagonal elements
    T3 = diag(a,-1)+diag(be)+diag(c,+1);
    T{3} = T3;              % Even cosine determinant matrix
    T4 = zeros(size(T3));
    T4(1,1) = 0 ;T4(1,2) = e ;T4(2,1) = 2*e;
    T4(2:end,2:end)=T3(1:end-1,1:end-1);
    T{4} = T4 ;             % Even sine determinant matrix
    % Calculate the Eigenvalues of the Determinant matrices
    for j = 1:4
        E = eig(T{j});
        EE{j}(:,i) = E;
    end
end
epsilon = epsilon';
E1 = EE{1}';
E2 = EE{2}';
E3 = EE{3}';
E4 = EE{4}';
E = [E1 E2 E3 E4];
figure;
plot(E,epsilon,'r','Linewidth',1);
xlabel('\delta');
ylabel('\epsilon');
hold on;
axis([-10,45,-20,20]);
PlotAxisAtOrigin;

function PlotAxisAtOrigin
%PlotAxisAtOrigin Plot 2D axes through the origin

% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./40;
Yoff=diff(get(gca,'YLim'))./40;

% DRAW AXIS LINEs
plot(get(gca,'XLim'),[0 0],'k');
plot([0 0],get(gca,'YLim'),'k');

% Plot new ticks
for i=1:length(X)
    plot([X(i) X(i)],[0 Yoff],'-k');
end
for i=1:length(Y)
   plot([Xoff, 0],[Y(i) Y(i)],'-k');
end

% ADD LABELS
t1 = text(X,zeros(size(X))-2.*Yoff,XL);
t2 = text(zeros(size(Y))-3.*Xoff,Y,YL);
set(t1,'FontName','TimesNewRoman','Fontsize',10) ;
set(t2,'FontName','TimesNewRoman','Fontsize',10) ;
% Add Axis Labels
l1 = text(2,max(Y),'a');
set(l1,'FontName','TimesNewRoman','Fontsize',10) ;
l2 = text(max(X),2.,'q');
set(l2,'FontName','TimesNewRoman','Fontsize',10) ;
% Add title
t = title('Stability chart of Mathieu Equation') ;
set(t,'FontName','TimesNewRoman','Fontsize',10) ;
box off;
% axis square;
axis off;
set(gcf,'color','w');
