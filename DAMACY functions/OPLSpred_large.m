function YHAT = OPLSpred_large(Xnew,P_o,W_o,W,Q)

% Johan Westerhuis
% Biosystems Data Analysis
% University of AMsterdam

% JAW: 20 juli 2007
%
% Edited by G.H. Tinnevelt 26 feb 2018

LV = size(W,2) -1 ; 
Enew = Xnew;

YHAT(:,1) = (inv(W(:,1)'*W(:,1))*Enew*W(:,1))*Q(:,1);


for lv = 1:LV
    tnew_o = inv(W_o(:,lv)'*W_o(:,lv))*Enew*W_o(:,lv);
    Enew = Enew - tnew_o*P_o(:,lv)';
    YHAT(:,lv+1) = (inv(W(:,lv+1)'*W(:,lv+1))*Enew*W(:,lv+1))*Q(:,lv+1);
end

