function F=Chen2010_ReturnFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)

F=-Inf;

% price of housing
p=(r+delta_r)/(1+r);
% wage
w=(1-alpha)*((r+delta_k)/alpha)^(alpha/(alpha-1));

% housing transactions costs
if abs(hprime-h)<1e-8
    % Housing does not change --> zero adjustment costs
    tau_hhprime=0;
else
    tau_hhprime=phi*h;
end

% earnings
earnings=w*kappaj*z;

%% Renter
if hprime==0
    % Budget constraint
    cspend=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime; % -hprime=0, so Chen (2010) omits it, but I leave it here
    % cspend=c+p*d (consumption goods plus housing services)
    % Analytically, we can derive the split of cspend into c and p*d as
    c=theta*cspend;
    d=(cspend-c)/p;
    if c>0
        % utility function
        %uinner=(c^theta)*(d^(1-theta)); 
        %F=(uinner^(1-sigma))/(1-sigma);
        F = f_util(c,d,theta,sigma);
    end
    % Impose borrowing constraints
    if aprime<0
        F=-Inf;
    end
%% Owner
elseif hprime>0
    % I don't get why Chen2010 sets Tr to only go to homeowners (unimportant comment)
    % Budget constraint
    c=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime;
    if c>0
        % utility function
        %uinner=(c^theta)*(hprime^(1-theta)); 
        %F=(uinner^(1-sigma))/(1-sigma);
        F = f_util(c,hprime,theta,sigma);
    end
    % Impose collateral constraints
    if aprime<-(1-gamma)*hprime % 'mortgage' cannot be greater than (1-gamma)*hprime (can be thought of as the downpayment requirement)
        F=-Inf;
    end
end %end if hprime=0

end %end function "Chen2010_ReturnFn"
