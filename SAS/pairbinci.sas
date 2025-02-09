********************************************************************************;
*
* Program Name   : PAIRBINCI.SAS
* Author         : Pete Laud, SSU, University of Sheffield
* Date Created   : 07 Feb 2025
* Version date   : 07 Feb 2025
* Repository     : Download latest version from https://github.com/PeteLaud/ratesci-sas
*
* Type           : Macro
* Description    : Computes two-sided confidence intervals (CI) for the difference
*                  between two paired binomial proportions 
*
* Macro            a,b,c,d = input values from 2x2 contingency table, 
* variables:       
*                  
*                  level = (2-sided) confidence level required, e.g. 0.95 for 95% CI
*                          NB this corresponds to a NI test at the (1-LEVEL)/2 
*                          significance level*;

* subroutine to evaluate score at a given value of x;
%macro fx(x); 
  %if contrast = "RD" %then %do;
    theta = &x.;
    AA = 2 * n;
    BB = -1*b - c + (2 * n - b + c) * theta;
    CC = -1*c * theta * (1 - theta);
    p21 = (SQRT(BB**2 - 4 * AA * CC) - BB)/(2 * AA);

	denterm = n * (2 * p21 + theta * (1 - theta));
	if denterm = 0 then fxi = 1E8 * sign(b - c - n * theta);
    else fxi = (b - c - n * theta) / SQRT(lambda * n * (2 * p21 + theta * (1 - theta)));

	* Skewness correction calculations;
    p12 = p21 + theta;
    if a = 0 then p11 = 0;
    else p11 = a / (a + d) * (1 - p12 - p21);
    p22 = 1 - p11 - p12 - p21;
    p2d = min(1, max(0, p21 + p11));
    p1d = p2d + theta;
  %end;
  %if contrast = "RR" %then %do;

  %end;

  V = max(0, (p1d * (1 - p1d) + p2d * (1 - p2d) -
    2 * (p11 * p22 - p12 * p21)) / n) * lambda;
  mu3 = (p1d * (1 - p1d) * (1 - 2 * p1d) +
      ((-1)**3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * (-1) * (p11 * (1 - p1d)**2 + p21 * p1d**2 - p1d * p2d * (1 - p1d)) +
      3 * ((-1)**2) * (p11 * (1 - p2d)**2 + p12 * p2d**2 - p1d * p2d * (1 - p2d))
      ) / (n**2);
  if abs(mu3) < 1E-10 then scterm = 0; * Avoids issues with e.g. x = c(1, 3, 0, 6);
  else scterm = mu3 / (6 * V**(3 / 2));
  As = scterm;
  Bs = 1;
  Cs = -(fxi + scterm);
  num = (-Bs + sqrt(max(0, Bs**2 - 4 * As * Cs)));
  if (skew = "FALSE" | scterm = 0) then fx = fxi - Z;
  else fx = num / (2 * As) - Z;

%mend fx;

* Subroutine to find lower or upper confidence limit;
%macro Tango (a, b, c, d, X0, X1, Z);
  *data tangodata;
  iteration = 1; * Initializing iteration counter;
  n = a + b + c + d;
  X0 = &X0;
  X1 = &X1;
  Z = &Z; 
  if bcf = "FALSE" then lambda = 1;
  else lambda = n / (n - 1);
  change = 1;
  do until (abs(change) < 1E-10);

  %fx(x0);
  fx0 = fx;
  %fx(x1);
  fx1 = fx;

  /*
*+---------------------------------------------------+
* Evaluate x0;
*+---------------------------------------------------+;
    theta = x0;
    AA = 2 * n;
    BB = -1*b - c + (2 * n - b + c) * theta;
    CC = -1*c * theta * (1 - theta);
    p21 = (SQRT(BB**2 - 4 * AA * CC) - BB)/(2 * AA);

	denterm = n * (2 * p21 + theta * (1 - theta));
	if denterm = 0 then fx0i = 1E8 * sign(b - c - n * theta);
    else fx0i = (b - c - n * theta) / SQRT(lambda * n * (2 * p21 + theta * (1 - theta)));

	* Skewness correction calculations;
    p12 = p21 + theta;
    if a = 0 then p11 = 0;
    else p11 = a / (a + d) * (1 - p12 - p21);
    p22 = 1 - p11 - p12 - p21;
    p2d = min(1, max(0, p21 + p11));
    p1d = p2d + theta;
    V = max(0, (p1d * (1 - p1d) + p2d * (1 - p2d) -
      2 * (p11 * p22 - p12 * p21)) / n) * lambda;
    mu3 = (p1d * (1 - p1d) * (1 - 2 * p1d) +
      ((-1)**3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * (-1) * (p11 * (1 - p1d)**2 + p21 * p1d**2 - p1d * p2d * (1 - p1d)) +
      3 * ((-1)**2) * (p11 * (1 - p2d)**2 + p12 * p2d**2 - p1d * p2d * (1 - p2d))
    ) / (n**2);
    if abs(mu3) < 1E-10 then scterm = 0; * Avoids issues with e.g. x = c(1, 3, 0, 6);
    else scterm = mu3 / (6 * V**(3 / 2));
    As = scterm;
    Bs = 1;
    Cs = -(fx0i + scterm);
    num = (-Bs + sqrt(max(0, Bs**2 - 4 * As * Cs)));
    if (skew = "FALSE" | scterm = 0) then fx0 = fx0i - Z;
    else fx0 = num / (2 * As) - Z;

*+---------------------------------------------------+
* Evaluate x1;
*+---------------------------------------------------+;
    theta = x1;
    AA = 2 * n;
    BB = -1*b - c + (2 * n - b + c) * theta;
    CC = -1*c * theta * (1 - theta);
    p21 = (SQRT(BB**2 - 4 * AA * CC) - BB) / (2 * AA);
	denterm = n * (2 * p21 + theta * (1 - theta));
	if denterm = 0 then fx1i = 1E8 * sign(b - c - n * theta);
    else fx1i = (b - c - n * theta) / SQRT(lambda * n * (2 * p21 + theta * (1 - theta)));

	* Skewness correction calculations;
    p12 = p21 + theta;
	if a = 0 then p11 = 0;
	else p11 = a / (a + d) * (1 - p12 - p21);
    p22 = 1 - p11 - p12 - p21;
    p2d = min(1, max(0, p21 + p11));
    p1d = p2d + theta;
    V = max(0, (p1d * (1 - p1d) + p2d * (1 - p2d) -
      2 * (p11 * p22 - p12 * p21)) / n) * lambda;
    mu3 = (p1d * (1 - p1d) * (1 - 2 * p1d) +
      ((-1)**3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * (-1) * (p11 * (1 - p1d)**2 + p21 * p1d**2 - p1d * p2d * (1 - p1d)) +
      3 * ((-1)**2) * (p11 * (1 - p2d)**2 + p12 * p2d**2 - p1d * p2d * (1 - p2d))
    ) / (n**2);
    if abs(mu3) < 1E-10 then scterm = 0; * Avoids issues with e.g. x = c(1, 3, 0, 6);
    else scterm = mu3 / (6 * V**(3 / 2));
    As = scterm;
    Bs = 1;
    Cs= -(fx1i + scterm);
    num = (-Bs + sqrt(max(0, Bs**2 - 4 * As * Cs)));
    if (skew = "FALSE" | scterm = 0) then fx1 = fx1i - Z;
    else fx1 = num / (2 * As) - Z;
*/

*  +----------------------------------------------------+;
    x2 = max(-1, min(1, x1 - (fx1 * (x1 - x0)/(fx1 - fx0)))); * update value of theta to evaluate;
    change = x2 - x1;
 *+-------------------------------------------------------------------------+
* Set x0 and x1 for next iteration of loop, and increment iteration counter;
*+-------------------------------------------------------------------------+;
*    output;
    x0 = x1;
    x1 = x2;
    iteration = iteration + 1;
  end;
*run;
%mend Tango;


%macro pairbinci(a, b, c, d, level, bcf=FALSE, skew=FALSE, contrast="RD");

  data tangodata(keep = a b c d level UCL LCL bcf skew);
    retain a b c d level LCL UCL;
	length skew $8. bcf $8.;
    a = &a;
    b = &b;
    c = &c;
    d = &d;
	contrast = "&contrast.";
    level = &level.;
    skew = "&skew.";
    bcf = "&bcf.";
    zval = probit(1 - (1 - &level)/2);
    %Tango(a, b, c, d, -1, 1, zval); 
    LCL = x2; 
    %Tango(a, b, c, d, -1, 1, -zval);
    UCL = x2;
    output;
  run;

*  proc print;
*  run;

%mend pairbinci; 

%pairbinci(1,1,7,12,level=0.95, bcf=TRUE, skew=TRUE);

/*
data tantest;
pi = constant('pi'); *3.14159265358979323846264338327950288;
 do i = 0 to 100;
	fi = round(tan(pi * (i/100 + 1) / 4), 1E-10);
	output;
 end;
run;
*/
