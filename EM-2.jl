
using StatsBase
using Distributions
function MJPBisect(ini,fin,m,E,Q,Qii,T1,T2,cond)
  tray=[T1 ini;T2 fin]
  t0=T2-T1
  t=t0/2
  t1=-t
  # Transition matrix of time lag t
  P=expm(t*Q)
  # We assign the probabilities of transitioning through each state
  prob=zeros(m)
  for i=E
    prob[i]=P[ini,i]*P[i,fin]
  end

  # We modify the probabilities if we are being conditioned to have at least 2 jumps
  if cond
    if ini==fin
      prob[ini]-=exp(-t0*Qii[ini])
    elseif Qii[ini]==Qii[fin]
      r=exp(t1*Qii[ini])*t*exp(t1*Qii[ini])*Q[ini,fin]
      prob[ini]-=r
      prob[fin]-=r
    else
      r=exp(t1*Qii[ini])
      r=r*(r-exp(t1*Qii[fin]))*Q[ini,fin]/(Qii[fin]-Qii[ini])
      prob[ini]-=r
      prob[fin]-=r
    end
    # We sample to see wich state we passed through
    mid=sample(E,weights(prob))

    # Once we have the state, we check to see in which manner we did so (number of jumps)
    if ini==mid
      if mid==fin
        r=exp(t1*Qii[ini])
        rc=P[ini,mid]-r
        (bef,aft)=sample([(0,2),(2,0),(2,2)],weights([r*rc,rc*r,rc*rc]))
      else
        r=exp(-t*Qii[ini])
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(0,2),(2,1),(2,2)],weights([r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    else
      if mid==fin
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        r1=exp(t1*Qii[mid])
        (bef,aft)=sample([(1,2),(2,0),(2,2)],weights([r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      else
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(1,1),(1,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    end
  else
    # We sample to see wich state we passed through
    mid=sample(E,weights(prob))

    # Once we have the state, we check to see in which manner we did so (number of jumps)
    if ini==mid
      if mid==fin
        r=exp(t1*Qii[ini])
        rc=P[ini,mid]-r
        (bef,aft)=sample([(0,0),(0,2),(2,0),(2,2)],weights([r*r,r*rc,rc*r,rc*rc]))
      else
        r=exp(t1*Qii[ini])
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(0,1),(0,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    else
      if mid==fin
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        r1=exp(t1*Qii[mid])
        (bef,aft)=sample([(1,0),(1,2),(2,0),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      else
        if Qii[ini]==Qii[mid]
          r=t*exp(t1*Qii[ini])*Q[ini,mid]
        else
          r=(exp(t1*Qii[ini])-exp(t1*Qii[mid]))*Q[ini,mid]/(Qii[mid]-Qii[ini])
        end
        if Qii[mid]==Qii[fin]
          r1=t*exp(t1*Qii[mid])*Q[mid,fin]
        else
          r1=(exp(t1*Qii[mid])-exp(t1*Qii[fin]))*Q[mid,fin]/(Qii[fin]-Qii[mid])
        end
        (bef,aft)=sample([(1,1),(1,2),(2,1),(2,2)],weights([r*r1,r*(P[mid,fin]-r1),(P[ini,mid]-r)*r1,(P[ini,mid]-r)*(P[mid,fin]-r1)]))
      end
    end
  end

  if bef==1
    if Qii[fin]>Qii[ini]
      tray=vcat(tray,[T1+rand(Truncated(Exponential(1/(Qii[fin]-Qii[ini])),0,t)) mid])
    elseif Qii[fin]<Qii[ini]
      tray=vcat(tray,[T1+t-rand(Truncated(Exponential(1/(Qii[ini]-Qii[fin])),0,t)) mid])
    else
      tray=vcat(tray,[T1+rand(1)*t mid])
    end
  elseif bef==2
    tray=vcat(MJPBisect(ini,mid,m,E,Q,Qii,T1,T1+t,true),tray)
  end
  if aft==1
    if Qii[fin]>Qii[ini]
      tray=vcat(tray,[T1+t+rand(Truncated(Exponential(1/(Qii[fin]-Qii[ini])),0,t)) fin])
    elseif Qii[fin]<Qii[ini]
      tray=vcat(tray,[T2-rand(Truncated(Exponential(1/(Qii[ini]-Qii[fin])),0,t)) fin])
    else
      tray=vcat(tray,[T1+t+rand(1)*t fin])
    end
  elseif aft==2
    tray=vcat(tray,MJPBisect(mid,fin,m,E,Q,Qii,T1+t,T2,true))
  end
  return(tray)
end
#####################################################
function MJPBridgeBisection(tray,Q)
  m=length(Q[1,1:end])
  n=length(tray[1:end,1])-1
  E=1:m
  Qii=-diag(Q)
  tray0=hcat([],[])
  for i = 1:n
    mat=MJPBisect(tray[i,2],tray[i+1,2],m,E,Q,Qii,tray[i,1],tray[i+1,1],false)
    mat=sortrows(mat)
    rng=hcat(mat[1,1],mat[1,2])
    for k = 2:length(mat[1:end,1])
      if mat[k,2]!=mat[k-1,2]
        rng=vcat(rng,mat[k,1:2])
      end
    end
    tray0=vcat(tray0,rng)
  end
  tray0=vcat(tray0,tray[end;1:2])
  return(tray0)
end
Q=[-4 2 1 1;
   2 -6 3 1;
   1 1 -4 2;
   1 3 2 -6]
tray0=[0.1 2;
      0.15 3;
      0.25 1;
      0.40 2;
      0.67 1;
      0.88 4;
      1.2 2]
B=MJPBridgeBisection(tray0,Q)
length(B[:,1])


function cuenta(trayecto)
  tra=trayecto
  tam=size(tra)[1]
  E=tra[1:end,2]
  n=size(unique(E),1)
  T=tra[1:end,1]
  R=zeros(n)
  N=zeros(n,n)
  for i in 1:(tam-1)
    delta=T[i+1]-T[i]
    R[E[i]]=R[E[i]]+delta
    N[E[i],E[i+1]]=N[E[i],E[i+1]]+1
  end
  return (R,N)
end

function EstimMCC(trayecto)
  tra= trayecto
  tam=size(tra)[1]
  E=tra[1:end,2]
  n=size(unique(E),1)
  T=tra[1:end,1]
  R=zeros(n)
  N=zeros(n,n)
  Q=zeros(n,n)
  P=zeros(n)
  for i in 1:(tam-1)
    delta=T[i+1]-T[i]
    R[E[i]]=R[E[i]]+delta
    N[E[i],E[i+1]]=N[E[i],E[i+1]]+1
  end
  Q=N./R
    for i in 1:n
        Q[i,i]=-sum(Q,2)[i]
    end
    P=N./sum(N,2)
  return (Q,P)
end

cuenta(tray0)
EstimMCC(tray0)

N=cuenta(tray0)[2]
Ri=cuenta(tray0)[1]
Es=cuenta(tray0)[2]





function EMQ(tray,NE,RE,ite)
    E=tray[1:end,2]
    n=size(unique(E),1) 
    N=zeros(n,n)
    R=zeros(n)  
    EN=zeros(n,n)
    ER=zeros(n)
    Q=zeros(n,n,ite)
    for k in 1:ite
        tra=MJPBridgeBisection(tray,Q[:,:,k])
        m=length(tra[:,1])
        N=cuenta(tra)[2]
        R=cuenta(tra)[1]
        for i in 1:n
            for j in 1:n
                if i!=j
                    EN[i,j]=mean(N[i,j])
                    ER[i]=mean(R[i])
                    Q[i,j,k]=EN[i,j]/ER[i]
                    Q[i,j,1]=mean(NE[i,j])/mean(RE[i])
                end
                
            end
            Q[i,i,k]=-sum(Q[:,:,k],2)[i]
            Q[i,i,1]=-sum(Q[:,:,1],2)[i]
        end

    end
    return (Q)
end

EMQ(tray0,Es,Ri,100)

datos=[1     4
2     1
3     4
4     3
5     3
6     2
7     3
8     1
9     2
10    3
11    4
12    1
13    1
14    4
15    3
16    2
17    3
18    4
19    1
20    3
21    3
22    1
23    3
24    1
25    4
26    2
27    1
28    1
29    1
30    1
31    1
32    2
33    4
34    1
35    4
36    3
37    3
38    3
39    3
40    3
41    4
42    1
43    3
44    1
45    3
46    2
47    1
48    2
49    1
50    3
51    3
52    3
53    1
54    4
55    2
56    2
57    3
58    1
59    2
60    1
61    3
62    3
63    2
64    2
65    3
66    2
67    3
68    4
69    3
70    3
71    2
72    3
73    3
74    3
75    2
76    3
77    4
78    4
79    3
80    3
81    2
82    1
83    2
84    3
85    2
86    1
87    3
88    1
89    2
90    3
91    2
92    1
93    1
94    2
95    1
96    3
97    4
98    4
99    2
100   1
101   2
102   1
103   2
104   1
105   3
106   1
107   4
108   3
109   1
110   4
111   3
112   3
113   4
114   2
115   3
116   3
117   3
118   4
119   2
120   2
121   3
122   4
123   2
124   4
125   1
126   3
127   4
128   1
129   1
130   2
131   1
132   1
133   4
134   2
135   2
136   1
137   3
138   1
139   1
140   3
141   4
142   3
143   4
144   2
145   1
146   3
147   2
148   4
149   3
150   1
151   1
152   1
153   1
154   1
155   1
156   1
157   1
158   1
159   3
160   2
161   2
162   4
163   2
164   4
165   2
166   3
167   1
168   3
169   2
170   2
171   1
172   1
173   4
174   3
175   3
176   4
177   3
178   4
179   1
180   2
181   2
182   2
183   1
184   1
185   1
186   3
187   2
188   2
189   1
190   1
191   3
192   1
193   2
194   1
195   4
196   2
197   4
198   3
199   1
200   3
201   2
202   3
203   2
204   1
205   3
206   4
207   2
208   1
209   2
210   2
211   2
212   1
213   3
214   3
215   1
216   4
217   1
218   3
219   3
220   3
221   2
222   1
223   1
224   4
225   3
226   2
227   3
228   2
229   4
230   3
231   3
232   1
233   1
234   2
235   4
236   4
237   1
238   3
239   1
240   3
241   3
242   3
243   3
244   1
245   4
246   2
247   3
248   2
249   3
250   3
251   2
252   2
253   1
254   4
255   3
256   2
257   3
258   2
259   3
260   3
261   3
262   1
263   1
264   3
265   1
266   4
267   4
268   4
269   2
270   1
271   1
272   1
273   2
274   3
275   4
276   3
277   3
278   3
279   4
280   1
281   1
282   3
283   3
284   2
285   4
286   3
287   3
288   2
289   1
290   3
291   1
292   3
293   4
294   3
295   4
296   3
297   1
298   3
299   4
300   3
301   3
302   2
303   2
304   1
305   3
306   3
307   4
308   4
309   2
310   4
311   2
312   2
313   4
314   3
315   3
316   3
317   3
318   4
319   4
320   1
321   4
322   2
323   2
324   2
325   1
326   3
327   3
328   3
329   4
330   3
331   2
332   1
333   4
334   2
335   3
336   1
337   2
338   1
339   2
340   2
341   3
342   2
343   4
344   2
345   3
346   1
347   3
348   2
349   1
350   2
351   1
352   4
353   3
354   3
355   2
356   3
357   3
358   2
359   2
360   3
361   3
362   3
363   4
364   4
365   2
366   1
367   1
368   4
369   2
370   4
371   1
372   2
373   1
374   3
375   1
376   4
377   4
378   1
379   2
380   2
381   1
382   3
383   2
384   2
385   3
386   1
387   1
388   4
389   1
390   2
391   3
392   3
393   2
394   1
395   2
396   1
397   3
398   3
399   3
400   1
401   2
402   3
403   1
404   1
405   1
406   3
407   2
408   4
409   1
410   4
411   3
412   1
413   1
414   2
415   2
416   1
417   2
418   2
419   2
420   1
421   3
422   3
423   2
424   1
425   1
426   4
427   4
428   2
429   2
430   4
431   4
432   2
433   1
434   3
435   3
436   3
437   3
438   1
439   1
440   1
441   2
442   4
443   1
444   2
445   2
446   2
447   2
448   4
449   1
450   4
451   3
452   4
453   2
454   4
455   3
456   2
457   3
458   4
459   3
460   3
461   3
462   4
463   4
464   3
465   1
466   2
467   2
468   3
469   2
470   4
471   3
472   2
473   2
474   1
475   3
476   1
477   3
478   2
479   2
480   2
481   1
482   2
483   1
484   3
485   1
486   4
487   2
488   3
489   2
490   1
491   2
492   3
493   1
494   4
495   3
496   3
497   2
498   2
499   3
500   1
501   3
502   3
503   3
504   1
505   3
506   1
507   1
508   1
509   2
510   1
511   4
512   2
513   3
514   2
515   3
516   1
517   3
518   2
519   4
520   3
521   2
522   1
523   1
524   3
525   3
526   1
527   2
528   3
529   2
530   2
531   4
532   1
533   3
534   4
535   1
536   4
537   2
538   2
539   1
540   4
541   2
542   3
543   4
544   3
545   3
546   1
547   2
548   1
549   2
550   4
551   1
552   2
553   4
554   3
555   4
556   3
557   1
558   3
559   2
560   2
561   1
562   3
563   2
564   2
565   2
566   3
567   2
568   3
569   4
570   2
571   1
572   3
573   2
574   1
575   2
576   1
577   3
578   4
579   1
580   1
581   1
582   3
583   3
584   4
585   4
586   1
587   1
588   2
589   4
590   3
591   3
592   2
593   1
594   1
595   2
596   3
597   4
598   1
599   3
600   4
601   2
602   2
603   1
604   1
605   1
606   4
607   3
608   4
609   1
610   3
611   1
612   1
613   1
614   2
615   1
616   3
617   3
618   2
619   3
620   2
621   3
622   3
623   2
624   3
625   3
626   3
627   1
628   3
629   3
630   1
631   3
632   3
633   4
634   4
635   4
636   3
637   3
638   4
639   1
640   3
641   3
642   3
643   1
644   4
645   2
646   1
647   1
648   4
649   2
650   3
651   2
652   3
653   3
654   1
655   4
656   4
657   1
658   3
659   1
660   3
661   1
662   4
663   4
664   4
665   3
666   4
667   4
668   3
669   3
670   2
671   3
672   4
673   2
674   3
675   1
676   4
677   4
678   2
679   3
680   2
681   2
682   2
683   1
684   4
685   1
686   2
687   2
688   4
689   1
690   3
691   4
692   4
693   1
694   4
695   3
696   3
697   1
698   2
699   3
700   1
701   2
702   4
703   1
704   3
705   1
706   2
707   3
708   3
709   2
710   2
711   3
712   2
713   4
714   2
715   3
716   3
717   2
718   4
719   1
720   3
721   2
722   4
723   1
724   1
725   4
726   4
727   3
728   1
729   1
730   1
731   2
732   1
733   4
734   1
735   3
736   4
737   3
738   1
739   1
740   3
741   3
742   4
743   1
744   3
745   2
746   2
747   2
748   1
749   1
750   3
751   4
752   1
753   4
754   2
755   3
756   4
757   2
758   1
759   4
760   2
761   4
762   4
763   3
764   2
765   1
766   3
767   4
768   3
769   2
770   1
771   4
772   1
773   1
774   4
775   4
776   3
777   2
778   3
779   2
780   4
781   1
782   4
783   2
784   1
785   3
786   3
787   1
788   1
789   3
790   1
791   1
792   4
793   3
794   4
795   3
796   2
797   4
798   3
799   1
800   3
801   3
802   1
803   2
804   3
805   3
806   1
807   2
808   1
809   3
810   1
811   4
812   1
813   1
814   4
815   1
816   1
817   1
818   2
819   4
820   1
821   1
822   3
823   1
824   1
825   2
826   2
827   2
828   4
829   4
830   3
831   4
832   4
833   3
834   1
835   4
836   4
837   1
838   2
839   1
840   2
841   4
842   1
843   3
844   3
845   3
846   1
847   3
848   1
849   4
850   2
851   2
852   4
853   3
854   3
855   3
856   4
857   4
858   1
859   2
860   3
861   3
862   3
863   1
864   3
865   2
866   4
867   1
868   3
869   2
870   4
871   3
872   3
873   2
874   3
875   4
876   4
877   3
878   2
879   2
880   3
881   3
882   1
883   3
884   4
885   3
886   2
887   4
888   4
889   1
890   2
891   1
892   3
893   1
894   4
895   3
896   4
897   3
898   4
899   1
900   2
901   4
902   3
903   1
904   4
905   2
906   3
907   3
908   1
909   1
910   2
911   2
912   3
913   1
914   2
915   1
916   2
917   3
918   4
919   3
920   3
921   2
922   3
923   1
924   1
925   1
926   4
927   2
928   1
929   4
930   3
931   1
932   2
933   1
934   3
935   2
936   3
937   1
938   3
939   4
940   3
941   2
942   4
943   3
944   3
945   3
946   2
947   3
948   3
949   1
950   4
951   3
952   3
953   3
954   1
955   3
956   2
957   1
958   3
959   4
960   3
961   3
962   3
963   2
964   1
965   4
966   2
967   2
968   1
969   4
970   2
971   3
972   2
973   1
974   1
975   2
976   3
977   3
978   3
979   3
980   4
981   2
982   3
983   4
984   4
985   2
986   1
987   4
988   2
989   4
990   4
991   4
992   4
993   2
994   1
995   2
996   3
997   4
998   1
999   3
1000  3
1001  2
]

EstimMCC(datos)
cuenta(datos)
R=cuenta(datos)[1]
Qin=EstimMCC(datos)[1]
Qe=cuenta(datos)[2]



EMQ(datos,Qe,R,1000)


