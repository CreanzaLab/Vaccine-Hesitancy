%% Under Review: The Cultural Evolution of Vaccine Hesitancy: Modeling the Interaction between Beliefs and Behaviors
%% medRxiv 2022.05.26.22275604; doi: https://doi.org/10.1101/2022.05.26.22275604 

% Code set to run changes in phenotipic frequences over the course of a simulation for varying initialization vectors
% "vars" is an initializing vector structured [alpha1, alpha2, b0, b1, b2, b3, C0, C1, C2, C3, x1, x2, x3, x4]
% vars used: FourCorn =[0    0    0.0100    0.5000    0.5000    0.9900    0.0100    0.5000    0.5000    0.9900    0.8100    0.1000    0.0700    0.0200
%                       0    0    0.0100    0.5000    0.5000    0.9900    0.0100    0.5000    0.5000    0.9900    0.1000    0.8100    0.0200    0.0700
%                       0    0    0.0100    0.5000    0.5000    0.9900    0.0100    0.5000    0.5000    0.9900    0.0700    0.0200    0.1000    0.8100
%                       0    0    0.0100    0.5000    0.5000    0.9900    0.0100    0.5000    0.5000    0.9900    0.0200    0.0700    0.8100    0.1000]

function xProp= PolySimsIterate_Vaccinesfx2_Phase1_Lines_Annotated(vars)
close all
format short

rows=size(vars,1);

%CC = [0.1, 0.5, 0.8]; % Used for multipanel confidence transmission

%% Maximum Cultural Selection Coefficient
max_select = 0.1; %[-0.1, 0, 0.1, 0.5] Vector used for multipanel maximum selection coefficient figures

%% Arrays used for parameter ranges when applicable
%b11 = 1:-0.04:0.00;
%s11 = -0.1:0.03:0.5;% s11 = [0.5:-0.03:-0.1] used if on the y-axis % 
%C11 = 0.01:0.05:1;


%% Commented loops below for single parameter simulations
%for p = 1:length(max_select) % Used for multipanel maximum selection coefficient figures
    
 %for t = 1:length(b11) % Vertical Axis
     
    %C1 = b11(t); % "C1 = C2" varied on y axis according to b11 array. C1 is set here and C2 is set equal to C1 below
      %B12 = b11(t);
      %alpha1 = b11(t);
      %max_select = s11(t);
    
    
%for q = 1:length(s11) %Horizontal Axis
     
     %max_select = s11(q); %Maximum Selection Coefficient varied according to s11 array
     %C1 = C11(q);
     %c1 = C11(q);
     %alpha2= C11(q);
     %B12 = C11(q);
 
 for j=1:rows
       
nGen = 300;
xProp = zeros(nGen,5);


%% Mating Assortment
 alpha1= vars(j,1);% A+ choosing parent
 alpha2= vars(j,2);% A- choosing parent


%Parent Vaccination State Influence (bm) (Part of the input vector)
b0 = vars(j,3);%V- x V-
b1 = vars(j,4);%V- x V+
b2 = vars(j,5);%V+ x V-
b3 = vars(j,6);%V+ x V+

%Belief Influence on Vaccination (cn)
c0 = 0.01;%A- x A-
c1 = 0.5; %A- x A+
c2 = 0.5; %A+ x A-
c3 = 0.99;%A+ x A+

%Probability of transmitting vaccine confidence (Cn). Unless varied, specified in the input vector
C0= vars(j,7);%A- x A-;
C1= vars(j,8);%A- x A+; 
C2 =vars(j,9);%A+ x A-;
C3= vars(j,10);%A+ x A+;


%Vaccination Probability
B03 = (c3*(1+b0))/2; 
B02 = (c2*(1+b0))/2;
B01 = (c1*(1+b0))/2; 
B00 = (c0*(1+b0))/2; 
B13 = (c3*(1+b1))/2;
B12 = (c2*(1+b1))/2;%* (Parameters with "*" set equal to one another when varying Vaccination Probabality)  
B11 = (c1*(1+b1))/2;%*
B10 = (c0*(1+b1))/2;
B23 = (c3*(1+b2))/2;
B22 = (c2*(1+b2))/2;%*
B21 = (c1*(1+b2))/2;%*
B20 = (c0*(1+b2))/2;
B33 = (c3*(1+b3))/2;
B32 = (c2*(1+b3))/2;
B31 = (c1*(1+b3))/2;
B30 = (c0*(1+b3))/2;


%Initial Phenotype Frequencies, given in input vector vars, put into xProp vector
xProp(1,1)= vars(j,11);%V+A+
xProp(1,2)= vars(j,12);%V+A-
xProp(1,3)= vars(j,13);%V-A+
xProp(1,4)= vars(j,14);%V-A-
    
xProp(1,5) = xProp(1,1)+xProp(1,2)+xProp(1,3)+xProp(1,4);

        
 for i = 2:nGen
              
    x1 = xProp(i-1,1); %V+A+
    x2 = xProp(i-1,2); %V+A-
    x3 = xProp(i-1,3); %V-A+
    x4 = xProp(i-1,4); %V-A-
    
 %% Vaccination Frequency (x1 + x2) from previous iteration
    V = zeros(1,length(nGen));
    V(i-1) = xProp(i-1,1) + xProp(i-1,2);
    
    %% Mating Frequencies based on A+/A- preference
    m11 = (x1^2)*(1 - alpha1) + (alpha1*(x1^2))/(x1 + x3);
    m12 = (x1*x2)*(1 - alpha1);
    m13 = (x1*x3)*(1 - alpha1) + (alpha1*(x1*x3))/(x1 + x3);
    m14 = (x1*x4)*(1 - alpha1);
    
    m21 = (x2*x1)*(1 - alpha2);
    m22 = (x2^2)*(1 - alpha2) + (alpha2*(x2^2))/(x2 + x4);
    m23 = (x2*x3)*(1 - alpha2);
    m24 = (x2*x4)*(1 - alpha2) + ((x2*x4)*alpha2)/(x2+x4);
       
    m31 = (x3*x1)*(1 - alpha1) + ((x3*x1)*alpha1)/(x3+x1);
    m32 = (x3*x2)*(1 - alpha1);
    m33 = (x3^2)*(1 - alpha1) + (alpha1*(x3^2))/(x3 + x1);
    m34 = (x3*x4)*(1 - alpha1);
    
    m41 = (x4*x1)*(1 - alpha2);
    m42 = (x4*x2)*(1 - alpha2) + ((x4*x2)*alpha2)/(x2 + x4);
    m43 = (x4*x3)*(1 - alpha2);
    m44 = (x4^2)*(1 - alpha2) + (alpha2*(x4^2))/(x2 + x4);
    
    %%%%%%%%%%%%%%%%%%%%%
   %% Vertical Transmission and Cultural Selection

   %% Vaccination frequency dependent cultural selection coefficient.
   
   %Constants of selection coefficient function
    k = 13;
    n = 2/(exp(k/2) - exp(-k/2));

    %Selection associated with V+ trait
    s1 = zeros(1,length(max_select));
    s1(i-1) = -(((0.3)/( 1 + exp(-k*(V(i-1)-0.9)))-n)-max_select);% max_select(p) used with multipanel figures
    % for k = 13, function reduces to format in doi: https://doi.org/10.1101/2022.05.26.22275604
    
   %% Proportion of each phenotype in population after vertical transmission
       
    %V+A+ (Vaccinated-Confident)
    wBarx1prime=(1+s1(i-1))*((m11*B33*C3+m12*B32*C2+m21*B31*C1+m13*B23*C3+m31*B13*C3+m14*B22*C2+m41*B11*C1+m22*B30*C0...
        +m23*B21*C1+m32*B12*C2+m24*B20*C0+m42*B10*C0+m33*B03*C3+m34*B02*C2+m43*B01*C1+m44*B00*C0));
    
    %V+A- (Vaccinated-Hesitant)
    wBarx2prime=(1+s1(i-1))*((m11*B33*(1-C3)+m12*B32*(1-C2)+m21*B31*(1-C1)+m13*B23*(1-C3)+m31*B13*(1-C3)+m14*B22*(1-C2)...
        +m41*B11*(1-C1)+m22*B30*(1-C0)+m23*B21*(1-C1)+m32*B12*(1-C2)+m24*B20*(1-C0)+m42*B10*(1-C0)+m33*B03*(1-C3)...
        +m34*B02*(1-C2)+m43*B01*(1-C1)+m44*B00*(1-C0)));
    
    %V-A+ (Unvaccinated-Confident)
    wBarx3prime= ((m11*C3*(1-B33)+m12*C2*(1-B32)+m21*C1*(1-B31)+m13*C3*(1-B23)+m31*C3*(1-B13)+m14*C2*(1-B22)+m41*C1*(1-B11)...
        +m22*C0*(1-B30)+m23*C1*(1-B21)+m32*C2*(1-B12)+m24*C0*(1-B20)+m42*C0*(1-B10)+m33*C3*(1-B03)+m34*C2*(1-B02)+m43*C1*(1-B01)...
        +m44*C0*(1-B00))); 
    
     %V-A- (Unvaccinated-Hesitant)
    wBarx4prime=((m11*(1-C3)*(1-B33)+m12*(1-C2)*(1-B32)+m21*(1-C1)*(1-B31)+m13*(1-C3)*(1-B23)+m31*(1-C3)*(1-B13)+m14*(1-C2)*(1-B22)...
        +m41*(1-C1)*(1-B11)+m22*(1-C0)*(1-B30)+m23*(1-C1)*(1-B21)+m32*(1-C2)*(1-B12)+m24*(1-C0)*(1-B20)+m42*(1-C0)*(1-B10)...
        +m33*(1-C3)*(1-B03)+m34*(1-C2)*(1-B02)+m43*(1-C1)*(1-B01)+m44*(1-C0)*(1-B00)));
        
          
     wBar1 = wBarx1prime + wBarx2prime + wBarx3prime + wBarx4prime;
     
     %Normalization: divide each phenotype frequency by the sum of all phenotype frequencies to ensure everything sums to 1
     
     sum2 = (wBarx1prime/wBar1) + (wBarx2prime/wBar1)+ (wBarx3prime/wBar1) + (wBarx4prime/wBar1);
    
     xProp(i,:) = [(wBarx1prime/wBar1) (wBarx2prime/wBar1) (wBarx3prime/wBar1) (wBarx4prime/wBar1) sum2]; 
     
    
 end % Generation end
   
  
%% Plot
 
  figure(1)
  subplot(2, 2, j)
 
    plot(1:nGen,xProp(:,1),'r','LineWidth',1)% V+A+
    hold on
    plot(1:nGen,xProp(:,2),'b','LineWidth',1)% V+A-
    plot(1:nGen,xProp(:,3),'g','LineWidth',1)% V-A+
    plot(1:nGen,xProp(:,4),'c','LineWidth',1)% V-A-
    plot(1:nGen,sum(xProp(:,1:2),2), 'k','LineWidth',1) %V+
    plot(1:nGen,sum(xProp(:,[1,3]),2), 'm','LineWidth',1) %A+
    xlim([1 nGen])
    ylim([0 1])
    title(['vars= ', num2str(xProp(1,1)),' , ',num2str(xProp(1,2)),' , ', num2str(xProp(1,3)),' , ',num2str(xProp(1,4))])
    ylabel('Phenotype Frequencies')
    legend('V^+A^+', 'V^+A^-', 'V^-A^+', 'V^-A^-', 'V^+', 'A^+')
    xlabel('Iterations')
   
  
 end % vars end

end % function end
%end
%end
