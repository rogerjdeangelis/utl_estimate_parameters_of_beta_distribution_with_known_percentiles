Estimate parameters of beta distribution with known percentiles

   WPS/PROC R or IML/R

   The SAS code(CDF Beta) used to check the WPS result worked in both WPS and SAS
   
   Recent addition and related solution using SAS/IML on end
   Rick Wicklin via listserv.uga.edu

https://tinyurl.com/y7e2byzu
https://github.com/rogerjdeangelis/utl_estimate_parameters_of_beta_distribution_with_known_percentiles

https://tinyurl.com/yb6aa7pe
https://stackoverflow.com/questions/51091293/how-do-i-estimate-parameters-of-beta-distribution-with-known-percentiles-in-r

Onyambu profile
https://stackoverflow.com/users/8380272/onyambu


INPUT
=====

  Given three  points on a Beta cumulative distribution function
  estimate the parameters of the resulting beta distribution


            Cumulative
      X     Prob < x

   0.264   .025
   0.511   .50
   0.759   .975

                BETA  PDF

              a-1            b-1
             x      ( 1 - x )
   f(x) =  ----------------------
                 B(a,b)

              _      _
             | (a)  | (b)
    B(a,b) = ------------
               _
              | (a + b)



           Beta probability density function a,b>1

PFY |
  3 +                           BBA
    |                          B  AB
    |                         B     B
    |                        B       A
    |                       A        A
    |                       A         A
    |                      A          AA
  2 +                      A           A
    |                     A             A
    |                     A             A
    |                    B               A
    |                   A                A
    |                   A                 A
    |                  A                  A
  1 +                  A                   A
    |                 A                    AA
    |                AA                     A
    |                A                       A
    |               B                        AA
    |             AB                          AB
    |           ABA                             BB
  0 +       ABBBA                                 BBB
    |
    ---+---------+---------+---------+---------+---------+--
      0.0       0.2       0.4       0.6       0.8       1.0

                                 X

  Cumulative beta probability density function a,b>1
  with given points

CDFY |                                        |
     |0.975                                   |
1.00 +---------------------------------------*******
     |                                     ****
     |                                    **  |
     |                                  **    |
     |      Given (x,y=prob < x)       **     |
0.75 +      3 points estimate beta    **      |
     |                               **       |
     |                              **        |
     |0.511                        **         |
     |----------------------------**          |
0.50 +                            *           |
     |                           *|           |
     |                          * |           |
     |                         *  |           |
     |                        **  |           |
0.25 +                       **   |           |
     |                      **    |           |
     |                    **      |           |
     |0.025              **       |           |
     |---------------|***         |           |
0.00 +       *********            |           |
     |               | 0.264      | 0.511     | 0.759
     ---+---------+---------+---------+---------+---------+--
       0.0       0.2       0.4       0.6       0.8       1.0

                          X

EXAMPLE OUTPUT
--------------

   a= 7.37502
   b= 7.07158

             6.4            6.1
            x      ( 1 - x )
   f(x) =  --------------------
               B(7.4,7.1)

              _      _
             | (7.4)  | (7.1)
    B(a,b) = ----------------
               _
              | (14.5)

    Check

           INPUT
        ==============       Estimates
         X       CUMX     Using fitted a,b

       0.264    0.025    0.02500
       0.511    0.500    0.50000
       0.759    0.975    0.97750


PROCESS  (Working R code)
==========================

* techinque is to minimize the difference between estimated ND known probabilties
  at the known x locations;


m=function(z){
   if(any(z<0))return(NA);
   x=z[1];
   y=z[2];

   a=c(qbeta(&cumx1,x,y),qbeta(&cumx2,x,y),qbeta(&cumx3,x,y))-c(&x1,&x2,&x3);
   b=c(pbeta(&x1,x,y),pbeta(&x2,x,y), pbeta(&x3,x,y))-c(&cumx1,&cumx2,&cumx3);

   abs(a)+abs(b)
 };

sol=optim(c(0.1,1),function(x)sum(abs(m(x))),m)$par;
sol;

* SAS/WPS check;
data check;
  set have;
    cdfest=CDF('BETA',x,7.37502,7.07158) ;
run;quit;



OUTPUT
======

WORK.WANT total obs=2

         SOL

 a     7.37502
 b     7.07158

WORK.CHECK total obs=3

    X       CUMX     CDFEST

  0.264    0.025    0.02500
  0.511    0.500    0.50000
  0.759    0.975    0.97750


*                _              _       _
 _ __ ___   __ _| | _____    __| | __ _| |_ __ _
| '_ ` _ \ / _` | |/ / _ \  / _` |/ _` | __/ _` |
| | | | | | (_| |   <  __/ | (_| | (_| | || (_| |
|_| |_| |_|\__,_|_|\_\___|  \__,_|\__,_|\__\__,_|

;

data have;
  x=0.264;  cumx=.025 ;call symputx('x1',x);call symputx('cumx1',cumx);output;
  x=0.511;  cumx=.50  ;call symputx('x2',x);call symputx('cumx2',cumx);output;
  x=0.759;  cumx=.975 ;call symputx('x3',x);call symputx('cumx3',cumx);output;
run;quit;

*          _       _   _
 ___  ___ | |_   _| |_(_) ___  _ __
/ __|/ _ \| | | | | __| |/ _ \| '_ \
\__ \ (_) | | |_| | |_| | (_) | | | |
|___/\___/|_|\__,_|\__|_|\___/|_| |_|

;


%utl_submit_wps64('
options set=R_HOME "C:/Program Files/R/R-3.3.2";
libname wrk  sas7bdat "%sysfunc(pathname(work))";
proc r;
submit;
source("C:/Program Files/R/R-3.3.2/etc/Rprofile.site", echo=T);
library(haven);
have<-read_sas("d:/sd1/have.sas7bdat");
head(have);
m=function(z){
   if(any(z<0))return(NA);
   x=z[1];
   y=z[2];
   a=c(qbeta(&cumx1,x,y),qbeta(&cumx2,x,y),qbeta(&cumx3,x,y))-c(&x1,&x2,&x3);
   b=c(pbeta(&x1,x,y),pbeta(&x2,x,y), pbeta(&x3,x,y))-c(&cumx1,&cumx2,&cumx3);
   abs(a)+abs(b)
 };
sol=optim(c(0.1,1),function(x)sum(abs(m(x))),m)$par;
sol;
endsubmit;
import r=sol   data=wrk.want;
run;quit;
data check;
  set wrk.have;
    cdfest=CDF("BETA",x,7.37502,7.07158) ;
run;quit;
proc print data=check;
run;quit;
');


Rick Wicklin via listserv.uga.edu
6:55 PM (15 hours ago)
 to SAS-L
This is sometimes called quantile-matching estimation (QME). Because
the quantiles involve the cumulative distribution function (CDF), the
equation does not usually have a closed-form solution, even when the
number of quantiles matches the number of parameters. In this problem,
you have more quantiles then parameters. The solution that Roger presents
(from 'onyambu') is an L1 (leas
t absolute value) solution, but many researchers prefer an L2 solution,
which is either least square or weighted least squares.

In SAS you can use PROC NLIN to produce either L2 solution. A discussion
of the general problem (m quantiles and k parameters, m > k) is discussed
in the blog post "Fit a distribution from quantiles," where you will also find working SAS code:
https://blogs.sas.com/content/iml/2018/03/07/fit-distribution-matching-quantile.html

Regards,
Rick Wicklin


