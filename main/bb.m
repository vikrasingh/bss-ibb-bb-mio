function [rParaOut,xbest,fbest,Xstar]=bb(D,d,A,b,c,lb,ub,diagInv,iPara,rPara,IotherPara,textfileIntermOut,IstopCondPara,targetfbest,toDebug)
% Date 08 Nov 2021
% Author: Vikram Singh
% The improved Branch and Bound algorithm from the reference
% "Somol,Pudil 2004 Fast Branch and Bound algorithm for optimal feature selection. IEEE Transactions on pattern analysis and machine intelligence".

% output: fbest is the criterion function value at the optimal point xbest, where xbest is the fitted ols coeff.
% Xstar is the best set of features
% rParaOut is a structure to save real output parameter for the algorithm
% rParaOut.nfuneval, rParaOut.nNodeVisit, rParaOut.nDC1, rParaOut.savecpu, rParaOut.stopCriteriaFlag, rParaOut.necConMaxVioInfQM

%  4 Sep 2023, adding the structure below to be passed to quadratic minimization package
%  S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = diag vector  of order nNonZeroEig x 1 of the D matrix 
%                  S.diagInv= 1/diag(A)  nx1 vector

   % Initialization
   tstart=cputime;
   k=1;
   X.candset=(1:D);
   psi=(1:D);
   r=D;
   fbest=inf; % initialization
   Xstar=0;
   stop=0; 
   Q=zeros(D-d+1,D,'uint16');  % Ordered set of features assigned to edges, Q=[Q_0; Q_1;. . . ; Q_k]
                                           % where each Q_k=[Q_(k,1) . . . Q_(k,qk)]
   J=1.2345*ones(D-d+1,D);  % vector of obj function value, J=[J_0; J_1;. . . ; J_k], where each J_k=[J_(k,1) . . . J_(k,qk)]
   q=zeros(1,D,'uint16'); % initialization
   gostep2=0; % flag to get to step2
   nfuneval=0; % to save the number of criterion function evaluations
   nNodeVisit=0;
   nDC1=0; % to save the no. of deletin conditions 1 happened, if fbest<Inf f(Node), then prune the tree, this is equivalent to out DC1
   stopCriteriaFlag=0; % flag to print in the output file
                       % = 0 (desired) means went over all the combinations
                       % = 6 (hard stop) means stop because of cputime limit
   activeAlgo=IotherPara(23);           % 19Dec23 find( QuadMinFunPara==1,1 );  % the active algo. to find Inf F                    
   timeLimit=abs(IstopCondPara(6));
   necConMaxVioInfQM=-inf; % initialize as vio. of the initial feasible sol.  ,to save max violation of the first order necessary cond. by quadratic min. solution for Inf F in the total run
   minSolTreeOption=1;  % whether to use minimum solution tree suggested in Yu B.,Yuan B. PR(1993) A more efficient Branch and Bound algorithm for feature selection.
                        % 1= yes, =0 no

   % main loop                     
   while stop==0
      
       if ((cputime-tstart)/60) > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
          [necConMaxVioInfQM,xbest,fbest]=olsfit(D,Xstar,A,b,c,lb,ub,diagInv,[],iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM); % extra call the get the regss. coeff for the best set of pred. Xstar
          savecpu=cputime-tstart;
          nfuneval=nfuneval+1;
          if toDebug>=1
              fprintf(textfileIntermOut,'BB stops because cputime limit exceeds %1.2f min.\n',savecpu/60);
              fprintf('BB stops because cputime limit exceeds %1.2f min.\n',savecpu/60);
          end
          stopCriteriaFlag=6;
          rParaOut.nfuneval=nfuneval;rParaOut.nNodeVisit=nNodeVisit;rParaOut.nDC1=nDC1;rParaOut.savecpu=savecpu;rParaOut.stopCriteriaFlag=stopCriteriaFlag;rParaOut.necConMaxVioInfQM=necConMaxVioInfQM;
          return;
       end

       if IstopCondPara(1)==1
           if fbest-targetfbest <=eps
              [necConMaxVioInfQM,xbest,fbest]=olsfit(D,Xstar,A,b,c,lb,ub,diagInv,[],iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM); % extra call the get the regss. coeff for the best set of pred. Xstar
              savecpu=cputime-tstart;
              nfuneval=nfuneval+1;
              if toDebug>=1
                  fprintf(textfileIntermOut,'BB stops because fbest - targetfbest <= 0');
                  fprintf( 'BB stops because fbest - targetfbest <= 0' );
              end
              stopCriteriaFlag=1;
              rParaOut.nfuneval=nfuneval;rParaOut.nNodeVisit=nNodeVisit;rParaOut.nDC1=nDC1;rParaOut.savecpu=savecpu;rParaOut.stopCriteriaFlag=stopCriteriaFlag;rParaOut.necConMaxVioInfQM=necConMaxVioInfQM;
              return;
           end
       end

      % Step 1.
      if gostep2==0
         [stop,necConMaxVioInfQM,nfuneval,psi,r,k,q,J,Q]=getCurrentNodeDescend(X,Q,J,q,k,r,psi,d,D,A,b,c,lb,ub,diagInv,iPara,rPara,activeAlgo,nfuneval,targetfbest,necConMaxVioInfQM,tstart,timeLimit,fbest);
      end
      if stop==1  % cputime limit has reached  
         [necConMaxVioInfQM,xbest,fbest]=olsfit(D,Xstar,A,b,c,lb,ub,diagInv,[],iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM); % extra call the get the regss. coeff for the best set of pred. Xstar
          savecpu=cputime-tstart;
          nfuneval=nfuneval+1;
          if toDebug>=1
              fprintf(textfileIntermOut,'BB stops because cputime limit exceeds %1.2f min.\n',savecpu/60);
              fprintf('BB stops because cputime limit exceeds %1.2f min.\n',savecpu/60);
          end
          stopCriteriaFlag=6;
          rParaOut.nfuneval=nfuneval;rParaOut.nNodeVisit=nNodeVisit;rParaOut.nDC1=nDC1;rParaOut.savecpu=savecpu;rParaOut.stopCriteriaFlag=stopCriteriaFlag;rParaOut.necConMaxVioInfQM=necConMaxVioInfQM;
          return;      
      end
      
      % Step 1* change for the minimum solution tree approach
      if minSolTreeOption==1
          if gostep2==0   
              if k~=(D-d)        % if not a leaf i.e. k~=D-d
                  Qkassign=[];
                  for idx=1:k
                      Qkassign=union( Qkassign,Q(idx,q(idx)) );
                  end
                  [necConMaxVioInfQM,~,Jc]=olsfit( D,1:D,A,b,c,lb,ub,diagInv,union( Qkassign,psi ) , iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM );
                  nfuneval=nfuneval+1;
                  if Jc<fbest
                      % update bound
                      Xstar=setdiff( 1:D,union( Qkassign,psi ) );
                      fbest=Jc;
                      if toDebug>=1
                          fprintf(textfileIntermOut,'fbest=%1.4f; #f evals=%d \n',fbest,nfuneval);
                          fprintf('fbest=%1.4f; #f evals=%d \n',fbest,nfuneval);
                      end
                  end
                  % Step 3.
                  [psi,r,q,Q]=cutOffTree(Q,q,k,r,psi);  % this pruning is due the min. solution tree approach, not due the DC1
                  nNodeVisit=nNodeVisit+1;
              end
          end
      end
      gostep2=0;
      
      % Step 2.
      if q(k)==0
         % Step 4.
         [stop,k,Q,X]=backtrack(X,Q,q,k); 
         if stop==1
            [necConMaxVioInfQM,xbest,fbest]=olsfit(D,Xstar,A,b,c,lb,ub,diagInv,[],iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM); % extra call the get the regss. coeff for the best set of pred. Xstar
            savecpu=cputime-tstart;
            nfuneval=nfuneval+1;
            rParaOut.nfuneval=nfuneval;rParaOut.nNodeVisit=nNodeVisit;rParaOut.nDC1=nDC1;rParaOut.savecpu=savecpu;rParaOut.stopCriteriaFlag=stopCriteriaFlag;rParaOut.necConMaxVioInfQM=necConMaxVioInfQM;
            return;
         end
         
         % Step 3.
         [psi,r,q,Q]=cutOffTree(Q,q,k,r,psi);
         gostep2=1;
      end
      
      if gostep2==0
          if J(k,q(k)) > fbest
              % Step 3.
              [psi,r,q,Q]=cutOffTree(Q,q,k,r,psi);
              nDC1=nDC1+1;
              gostep2=1;
          else
              X(k+1).candset=setdiff(X(k).candset,Q(k,q(k)),'stable');
          end
      end
      
      if gostep2==0
          if k==(D-d)
              % Step 5.
              [Xstar,fbest]=updateBound(X,J,q,k);
              % Step 3.
              [psi,r,q,Q]=cutOffTree(Q,q,k,r,psi);
              gostep2=1;
          else
              k=k+1;
              nNodeVisit=nNodeVisit+1;
          end
      end
      
   end %===== end while stop=0
    
    
end %============ end BB(...)

%==================================================================================================================================================================================
function [stop,necConMaxVioInfQM,nfuneval,psi,r,k,q,J,Q]=getCurrentNodeDescend(X,Q,J,q,k,r,psi,d,D,A,b,c,lb,ub,diagInv,iPara,rPara,activeAlgo,nfuneval,targetfbest,necConMaxVioInfQM,tstart,timeLimit,fbest)
   % Step 1
   stop=0; % initialization
   q(k)=r-(D-d-k);
   % if using in-level node ordering
   Jk=zeros(1,length(psi)); % to store the value of the criterion function after removing every single feature
   for i=1:length(psi)
      [necConMaxVioInfQM,~,Jk(i)]=olsfit(D,X(k).candset,A,b,c,lb,ub,diagInv,psi(i),iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM);
      nfuneval=nfuneval+1;
      if ((cputime-tstart)/60) > timeLimit  % if the cputime becomes greater than timeLimit=|IstopCondPara(6)| min, then stop
         if fbest<inf
            stop=1;return
         end
      end
   end
   [Jk,Qk]=sort(Jk,'descend'); % sort from max to min for min J case

   psi=psi(Qk); % sort psi
   Q(k,1:q(k))=psi(1:q(k));  % to store index of qk no. of features after sorting in decreasing order
   J(k,1:q(k))=Jk(1:q(k));  % to store the value of the criterion function of qk no. features
   
   psi=setdiff(psi,Q(k,1:q(k)),'stable'); % exclude features from psi to avoid duplications
   r=r-q(k);
   
end

%==================================================================================================================================================================================

function [psi,r,q,Q]=cutOffTree(Q,q,k,r,psi)
   % Step 3
   
   psi=union(psi,Q(k,q(k)));
   r=r+1;
   Q(k,q(k))=0; % Q_k=Q_k\Q(k,q_k)
   q(k)=q(k)-1;

end

%==================================================================================================================================================================================

function [stop,k,Q,X]=backtrack(X,Q,q,k)
   % Step 4
   stop=0;
   k=k-1;
   if k==0
      stop=1;return
   else
      X(k).candset=union(X(k+1).candset,Q(k,q(k)));  
   end

end

%==================================================================================================================================================================================

function [Xstar,Fstar]=updateBound(X,J,q,k)
   % Step 5
   
   Fstar=J(k,q(k));
   Xstar=X(k+1).candset;
   
end
    
%==================================================================================================================================================================================

function [necConMaxVioInfQM,x_tilde,f_tilde]=olsfit(n,availInd,A,b,c,lb,ub,diagInv,dropInd,iPara,rPara,activeAlgo,targetfbest,necConMaxVioInfQM)

   isTF=0; % means do not use targetfbest for quadratic min. algo.
   keepIdx=setdiff(availInd,dropInd); % indices to keep to get the reduced dim. problem 
   nr=length(keepIdx); % reduced dim.
   lbr=lb(keepIdx);ubr=ub(keepIdx); % bounds of the box in the reduced dim.
   x0=0.5*(lbr+ubr); % starting point, just take the midpoint of the box
   
   % 5Sep23
   [Q,diagD]=eig(A(keepIdx,keepIdx),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
   Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
   Sstruct.Q=Q(:, (nr-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
   Sstruct.D=diagD( (nr-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
   Sstruct.diagInv=diagInv(keepIdx);

   [epsMax,xstar,f_tilde]=quad_min_algo_pool(nr,A(keepIdx,keepIdx),b(keepIdx),c,lbr,ubr,x0,Sstruct,iPara,rPara,activeAlgo,1,isTF,targetfbest);   % isSoftStop=1 hard coded
   x_tilde=zeros(n,1); x_tilde(keepIdx)=xstar; % get solution x in the original space
   if necConMaxVioInfQM<epsMax, necConMaxVioInfQM=epsMax;end

end

%==================================================================================================================================================================================

