function [dS] = yeastglycolysisNM(t,S)

S_1 = S(1); % Glucose
S_2 = S(2); % phosphate pool
S_3 = S(3); % 1,3-bisphosophoglycerate
S_4 = S(4); % Cytosolic pyruvate/ acetaldehyde pool
S_5 = S(5); % extracellular concentration S4 ^
A_3 = S(6); % ATP
N_2 = S(7); % NADH


% total ADP+ATP
A_tot = 4;
A_2 = A_tot -A_3;

% total NADH + NAD+
N_tot = 1;
N_1 = N_tot -N_2;

% incoming flux of glucose
J_G = 2.5; %mM/min
% J_G = 0.2;
% Hexokinase phosophoglucoisomerase and phosphofructokinase cooperatively
% inhibited with ATP 
k_1 = 100;  %mM/min max rxn rate
K_I = 0.52; % mM inhibition constant (related to half max)
q = 4;      % q is cooperativity coefficient   
v_1 = k_1*S_1.*A_3./(1+(A_3/K_I).^q); % mM/min rxn velocity
%v_1 = k_1*S_1*A_3;

% Glyceraldehydes-3-phosphate dehydrogenase linear in phosphate and ADP
k_2 = 6.0; %mM/min rxn constant
v_2 = k_2*S_2.*N_1; % mM^3/min <WTF!!!!!!! [k_4] must be 1/(min* mM^2)

% Phosphoglycerate kinase, phosophglycerate mutase, enolase, and pyruvate
% kinase
k_3 = 16.0; % mM/min rxn rate
v_3 = k_3*S_3.*A_2; % mM^3/min <WTF!!!!!!! [k_4] must be 1/(min* mM^2)

% Alcohol dehydrogenase
k_4 = 100; %mM/min
v_4 = k_4*S_4.*N_2; % mM^3/min <WTF!!!!!!! [k_4] must be 1/(min* mM^2)

% Nonglycolytic ATP consumption
k_5 = 1.28; % 1/min from supplement
% k_5 = 18; % infered from main text
v_5 = k_5*A_3; %mM/min 

% Formation of glycerol from triose phosphates
k_6 = 12.0; %mM/min supplement value
% k_6 = 6.0; % mM/min text value Table 3.
v_6 = k_6*S_2.*N_2; % mM^3/min <WTF!!!!!!! [k_6] must be 1/(min* mM^2)

% Degredation of pyruvate and acetaldehyde in the extracellular space
% k = 1.8; % 1/min in supplement
 k = 18; % infered from text 
v_7 = k*S_5;

% Carbon sink term to pyruvate pool acounting for carbon loss to cellular
% synthetic processes (Over specified model only)

% Membrane transport of pyruvate and acetaldehyde into extra cellular space
ASPV = 13.0; % 1/min
J_P = ASPV*(S_4 - S_5); % mM/min

% chemical species derivatives
phi = 0.10;


dS = zeros(7,1); % column vector
dS(1) = J_G*ones(size(v_1)) -v_1; % GLucose S1, S6
dS(2) = 2*v_1 - v_2 - v_6; % phosphate pool S1, S6, S2, S7
dS(3) = v_2 -v_3; % 1,3-bisphophoglycerate, S2, S6, S3, S7
dS(4) = v_3 -v_4 -J_P; % cytosolic pyruvate and acetaldehyde pool
dS(5) = phi*(J_P - v_7); % extracellular concentration of S_4 see above^
dS(6) = -2*v_1 + 2*v_3 - v_5; %ATP (A3) S1, S6, S3
dS(7) = v_2 -v_4 - v_6; %NADH (N2)



