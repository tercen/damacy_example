function [W, Q, W_o, P_o] = OPLS_large(X,y,LV)

% Johan Westerhuis
% Biosystems Data Analysis
% University of AMsterdam

% LV orthogonal components are calculated after which 1 predictive
% component is calculated

% JAW: 20 juli 2007

% Edited by G.H. Tinnevelt



P_o = zeros(size(X,2),LV);
W_o = zeros(size(X,2),LV);
W = zeros(size(X,2),LV+1);
Q = zeros(1,LV+1);
E = X;


for lv = 1:LV
    w = ((y'*y)\y'*E)';  
    w = w  / sqrt(w'*w);
    t = (w'*w)\E*w; 
    q = (t'*t)\t'*y;
    % c = (inv(t'*t)*t'*y)'; 
    % u = (inv*c'*c)*y*c;
    % These are two strange lines. For single y c is a scaler (so why
    % transposing it). Furthermore, c and u are not used anymore
    p = ((t'*t)\t'*E)'; 
    w_o = p - ((w'*w)\w'*p)*w; 
    w_o = w_o / sqrt(w_o'*w_o);
    t_o = (w_o'*w_o)\E*w_o; 
    p_o = ((t_o'*t_o)\t_o'*E)';

    E = E - t_o*p_o';

    P_o(:,lv) = p_o;
    W_o(:,lv) = w_o;
    W(:,lv) = w;
    Q(:,lv) = q;
end

w = ((y'*y)\y'*E)'; 
w = w / sqrt(w'*w);
W(:,LV+1) = w;
q = (t'*t)\t'*y;
Q(:,LV+1) = q;