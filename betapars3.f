* Use fit results from reweightfitg-0-0-Alt1.out to GP Run5
* Use reframe.f to change variables back to (alpha,beta)
*
* peak, body probabilities
      data pnorm/0.19189d0,0.29734d0/
      
* body beta parameters (alpha, beta, alpha-1, beta-1)
* Input (mean,rms) values    2.2823011253346112E-002   3.8547719901657680E-002
* beta   0.31972522119659613        13.689172099472287      -0.68027477880340381        12.689172099472287  
      data betabody/0.31973d0,13.689d0,-0.68027d0,12.689d0/   
*      data betabody/0.36074d0,14.707d0,-0.63926d0,13.707d0/

* arms beta parameters (alpha, beta, alpha-1, beta-1)
* Input (mean,rms) values    1.5837433426272810E-002   2.9153096463069059E-002
* beta   0.27460968197536101        17.064669643411381      -0.72539031802463905        16.064669643411381
      data betaarms/0.27461d0,17.065d0,-0.72539d0,16.065d0/
*      data betaarms/0.34137d0,18.161d0,-0.65863d0,17.161d0/
