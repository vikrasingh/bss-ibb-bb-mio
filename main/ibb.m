function   [numOfBoxDel,outputPara,ixtilde,ifxtilde,ytilde,Ytilde,xbest,fbest]=ibb(nPts,n,yArray,xMatrix,trueb,targetfbest,Anorm,bnorm,cnorm,X,lbX,ubX,diagInv,tmax,normxRelaxedOpt,prexbest,...
        normFactor,IstopCondPara,IotherPara,iPara,rPara,delCondPara,saveIntermOutput,numOfBoxDel,InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr,textfileName,toDebug,level8dirpath)
% 4 March 24, For Both AG71 and AG72 
% Implemented an acceleration strategy to immediately partition the current box if fbest improves at an iteration 
% 3 Jan24, For BrFS selection uses dlnode, for BFS,DFS,Max 1 and Max 0 selection uses HeapList
% 27 Nov23, introduced IotherPara(9)=-8m partition strategy 
% 19 Nov23, saving xbest, fbest in a matlab data file to be used as warm start for the next run

% 4 Sep 2023, adding the structure below to be passed to quadratic minimization package
%  S is a structure S.Q = Q matrix from spectral decomposition such that A=QDQ'  where Q is an orthonormal matrix and D is diagonal with eig. value in the decreasing order 
%                  S.nNonZeroEig = no. of non zero eigenvalues  
%                  S.D = diag vector  of order nNonZeroEig x 1 of the D matrix 
%                  S.diagInv= 1/diag(A)  nx1 vector

        % Initialize ouput parameters
        tStart=cputime;
        A_normdiv2=0.5*Anorm;  % as we need 0.5*A repeatedly to evaluate fx, we can define it once here to save cputime
        stop=0;          % will be used to check if the stopping condition holds or not.
        rankX=rank(xMatrix);
        quadMinAlg=IotherPara(23);  % quad. min. alg. to be used for data fitting
        [~,decWidIdx]=sort(ubX-lbX,'descend'); % indices of the width of a box from large to small,  to be used in selection of direction to partition the box
        
        % the following lines of code is to save the intermediate results in a text file
        if toDebug>=1
            intermFbestIter=saveIntermOutput.FbestIter; % initialization, using it to save the intermediate results at the given CPU time limits saveIntermOutput. Only when IstopCondPara(6)= -ive
            intermMaxIter=saveIntermOutput.MaxIter;
            intermCPU=saveIntermOutput.CPU;
            intermCPUCnt=1;intermFbestIterCnt=1;intermMaxIterCnt=1;
            nintermFbestIter=length(intermFbestIter);nintermMaxIter=length(intermMaxIter);nintermCPU=length(intermCPU); 
        end

        % set up direction to cut and no. of cuts for the child boxes
        ndigits=length( num2str(abs(IotherPara(9))) );
        cutDir=ceil(IotherPara(9)/(10^(ndigits-1)));         % first digit of the parameter will give the cut direction,i.e. which comp. to cut
        nCuts=abs( IotherPara(9)-cutDir*(10^(ndigits-1)) );  % no. of cuts
        if nCuts==0
            nCuts=1;
        end
        
        % set up box to be used for refining a feasible point and no. of iter. for refinement
        box_for_refine=floor( IotherPara(21)/(10^(length(num2str(IotherPara(21)))-1)) );         % first digit of IotherPara(21) will tell use which box to use for refinement
        niter_for_refine=IotherPara(21)-box_for_refine*(10^(length(num2str(IotherPara(21)))-1)); % remaining digits will tell us the no. of iterations to use for refinement
        modify_iPara=iPara;modify_iPara(2)=niter_for_refine;
        iParaGA=[];rParaGA=[];

        % set up selection critera of a box from the list
        firstSelOpt=floor( IotherPara(5)/(10^(length(num2str(IotherPara(5)))-1)) );    % first digit of IotherPara(5) defines primay selection criteria
        secndSelOpt=IotherPara(5)-firstSelOpt*(10^(length(num2str(IotherPara(5)))-1)); % remaining digits defines secondary selection criteria
          
        % save the list parameter in the matlab data file to continue with the algo. later on if we want to
        readatafile=fullfile(level8dirpath,sprintf('p%dd%de%ds%dt%d.mat',InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr));   % path to the matlab data file to store the algo. output para for resuming the algo.
        savedatafile=fullfile(level8dirpath,'saveListData',sprintf('p%dd%de%ds%dt%d.mat',InpFilCtr,pDimCtr,egCtr,setCtr,paraCtr));   % path to the matlab data file to store the algo. output para for resuming the algo.
        isFileExists=exist(readatafile,'file');          % if isFileExists=2 then yes, otherwise not

        % save cputime for intermediate steps
        cpuArray=zeros(1,10); 
        %{
         (1) cputime for branch function
         (2) cputime for expanding the list
         (3) cputime for selection and deletion of the box
         (4) cputime for selecting support for feasible pt subroutine
         (5) cputime for f evaluation
         (6) cputime for feasibility sampling
         (7) cputime for idx2flag setdiff
         (8) cputime for bisecting the parent box
         (9) cputime to process leaf node boxes using DC2 and DC5 
         (10) cputime for last fbest update 
        %}
        
        % Initialize ouput parameters
        nIter=0;         % initialize the num of iterations, which will keep on updating with each iteration.
        nBoxDeleted=0;   % will keep an account of how many boxes have been deleted from the list L
        nFeval=0;        % will give us the number of inclusion function evaluation for the run
        nfEval=0;        % will save the number of fx evaluation 
        nRefine=0;       % no. of refinements done for a feasible point
        necConMaxVioRefQM=-inf; % initialize, max. violation of the necessary cond. box or unbox for quadratic min. while doing refinement
        necConMaxVioInfQM=-inf; % initialize, max. violation of the necessary cond. box or unbox for quadratic min. while finding lb F
        cpuIncFunEval=0;
        cpuForRefineFesPt=0; % to save the cputime spent on refinement of the feasible point
        outputPara=(-1)*ones(18,1);  % will store the output parameter from the interval algo.
        outputPara(11)=0;  % initialization
        %{
 outputPara(1)= numOfBoxDeleted
            (2)=numOfBoxInlistL
            (3)=cputime for selecting and deleting nodes from the list
            (4)=numOfInc_FunEval
            (5)=numOfIter
            (6)=cputime for sampling
            (7)= which stopping criteria used
            (8)= saves the constraint violation for Xstar, but in setupForIntvalAlgo subroutine.
            (9)= cpu of inclusion function call
            (10)= number of zeros less than a tolerance in the solution Xstar
            (11)= last iteration at which fbest got updated
            (12)= num of small f calls
            (13)= cpu for sampling
            (14)= no. of refinements of the feasible point
            (15)= total no. of quadratic min. calls, i.e. sum of nFeval+nRefine, evaluated only at the end
            (16)=necConMaxVioRefQM   max. violation of the necessary cond. box or unbox for quadratic min. while doing refinement
            (17)=necConMaxVioInfQM   max. violation of the necessary cond. box or unbox for quadratic min. while finding lb F
            (18)=cputime in seconds for the algorithm
        %}
        
        % find an initial feasible point 
        if isFileExists==2  % we want to use the xbest,fbest from the last run as a warm start
           loadout=load(readatafile); % load the parameters to a structure
           ixtilde=loadout.xbest;ifxtilde=loadout.fbest;
        else 
            ifxtilde=inf;
            if ~isempty(prexbest) %------ if there is an initial solution available
                supprexstar=find(prexbest); % find the non zero entries
                nsupprexstar=length(supprexstar); % no. of entries in the supp
                if nsupprexstar==tmax  % prexstar is already feasible
                   ixtilde=prexbest;ifxtilde=feval(ixtilde(supprexstar),A_normdiv2(supprexstar,supprexstar),bnorm(supprexstar),cnorm); 
                else % if prexstar is not feasible
                    if nsupprexstar>tmax  % truncate the supp(prexstar)-tmax no. of components with smallest magnitude to get a feasible point
                       [~,idxKLargeEntLocal]=maxk(abs(prexbest(supprexstar)),tmax);
                       supprexstar=sort(supprexstar(idxKLargeEntLocal)); 
                       ixtilde=zeros(n,1); ixtilde(supprexstar)=prexbest(supprexstar); 
                    else % nsupprexstar<tmax , constuct the remaining support using xRelaxedopt
                       ctr=0; ictr=1;
                       while ctr<(tmax-nsupprexstar)
                          if ~ismember(decWidIdx(ictr),supprexstar) 
                             ctr=ctr+1; 
                             supprexstar=custom_union(supprexstar,decWidIdx(ictr));
                          end 
                          ictr=ictr+1;
                       end  
                       ixtilde=zeros(n,1); ixtilde(supprexstar)=prexbest(supprexstar);
                    end 
                    tempdegIntvls=1:n;tempdegIntvls(supprexstar)=[];
                    savecpu1=cputime;        
                    [necConMaxVioInfQM,ixtilde,ifxtilde]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,ixtilde,diagInv,targetfbest,tempdegIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                    cpuIncFunEval=cpuIncFunEval+(cputime-savecpu1);
                end % 
            end % endif ------------------   
           
            % initial feasible point
            [necConMaxVioRefQM,nRefine,cpuForRefineFesPt,cpuArray,supp,ixtilde2]=getFeasiblePt(n,[],1:n,lbX,ubX,normxRelaxedOpt,targetfbest,tmax,IotherPara,IstopCondPara,Anorm,bnorm,cnorm,diagInv,...
                                         rPara,cpuArray,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara);
            tempdegIntvls=1:n;tempdegIntvls(supp)=[];
            savecpu1=cputime;        
            [necConMaxVioInfQM,ixtilde2,ifxtilde2]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,ixtilde2,diagInv,targetfbest,tempdegIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
            cpuIncFunEval=cpuIncFunEval+(cputime-savecpu1);
            if ifxtilde2<ifxtilde
               ixtilde=ixtilde2;ifxtilde=ifxtilde2; 
            end
        end  % end if file exists

        % initialize xbest, fbest
        xbest=ixtilde;fbest=ifxtilde;
        
        % if using debugging mode
        if toDebug>=1
            disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
            fprintf('      ;Iteration; fbest; numOfBoxesInList; cputime in min.\n'); % to output in the command window
            fprintf('Initial= 0; %f; 1; 0\n',ifxtilde); 
            fprintf(textfileName.interm_out,'      ;Iteration; fbest; numOfBoxesInList; cputime in min.\n'); % to save in the text file
            fprintf(textfileName.interm_out,'Initial= 0; %f; 1; 0\n',ifxtilde);
        end

        % initial lbf(X)
        savecpu1=cputime; 
        [necConMaxVioInfQM,xlbfX,lbfX]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,ixtilde,diagInv,targetfbest,[],iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
        cpuIncFunEval=cpuIncFunEval+(cputime-savecpu1);
        nFeval=nFeval+1;

        % initialize the list
        nb=1;          % number of boxes will get updated throughout the run
        deleteElt=1;   % select the first box to start bisecting
        headnode=[];lastnode=[];
        if firstSelOpt==1 % (BrFS) box selection based on age of the box criteria using Linked List
           headnode=dlnode(X,lbfX,xlbfX,0,0,[],[]);   % assign the first node as head node
           lastnode=headnode;     % assign the last node as the head node as well for initialization
           heap=[];
           isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
        else  % for other selection options use HeapList
            nInitialBox=100;
            if firstSelOpt==4  % (BFS) box selection based on the lowest lb F(V)
               heap = minheap(nInitialBox, @(a,b) a.LBFBOX - b.LBFBOX, createbox );
               isSelBFS=1;    % 1 if the selected box is using min lbF, =0 if any other selection option
            elseif firstSelOpt==5  % (DFS) select box with max #2 flags
               heap = minheap(nInitialBox, @(a,b) b.N2FLAG - a.N2FLAG, createbox ); 
               isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
            elseif firstSelOpt==6  % select box with max #1 flags
               heap = minheap(nInitialBox, @(a,b) b.N1FLAG - a.N1FLAG, createbox );
               isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
            elseif firstSelOpt==7  % select box with max #0 flags
               heap = minheap(nInitialBox, @(a,b) b.N0FLAG - a.N0FLAG, createbox );
               isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
            end  % endif IotherPara(5)=1
            heap.add(createbox(X,lbfX,xlbfX,0,0,n,[],[])); % add the first box to the heap
        end

        ytilde=lbfX;Ytilde=X; % save output parameters
        lastfbest=fbest;saveIter=0;lastfbestcpu=cputime;isStop_dfbestCPUHolds=0;  % to use termination criteria based on improvement in fbest.  
        preSwitchIter=0;     % to control switching of selection criteria
        if n<=100, switchAfterIter=1000;  % no. of iter. after which we switch the selection criteria if there is no improvement in fbest
        elseif n<=300, switchAfterIter=200;
        elseif n<=500, switchAfterIter=100;
        elseif n<=1000, switchAfterIter=50;
        else, switchAfterIter=20;
        end
        switchSelOpt=0;   % =0 or 1 ,flag to control switching option for selection,   
        nBestBoxForSwitch=3; % compare the best given no. of boxes along with
        nRandBoxForSwitch=10;  % the given no. of randomly selected boxes to select a box with secondary selection criteria
        if toDebug==2,fprintf(textfileName.box_flag,'tmax=%d \n',tmax);end

        % main loop for the algo.
        while stop==0  % loop to use recursion until the stopping criteria does not get satisfied
            
            nIter=nIter+1;    % update num of iterations after each recursion.

            % Check stopping conditions
            if nb==0    % if the list become empty, i.e. all the combinations has been processed
                if toDebug>=1
                   fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,(cputime-tStart)/60);disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                   disp('Algo. 7 for BSS has stopped. All the combinations has been processed.');
                end
                outputPara(7)=0;
                break;  % stop while loop
            end

            ctime=(cputime-tStart)/60;    
            if ctime >= abs(IstopCondPara(6)) % if CPU time exceeded the hard limit provided
                if toDebug>=1
                    fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);fprintf(textfileName.interm_out,'Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);
                    disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                    disp('Termination criteria cpu time limit has reached holds');
                end
                outputPara(7)=6;
                break ; % exit while loop
            end

            if IstopCondPara(5)~=0   % IstopCondPara(5)=1 check this condition, IstopCondPara(5)=0, skip it
                if nIter==abs(IstopCondPara(5)) % if no. iter. reaches the max. iter. limit
                    ctime=(cputime-tStart)/60;
                    if toDebug>=1
                        fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);fprintf(textfileName.interm_out,'Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);
                        disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                        disp('Max. no. of iterations reached.');
                    end
                    outputPara(7)=5;
                    break  % exit while loop
                end
            end

            if IotherPara(3)~=0  % IotherPara(3)=1 check this condition, IotherPara(3)=0, skip it
               if isStop_dfbestCPUHolds==1 || ((cputime-lastfbestcpu)/60)>IstopCondPara(11) % if fbest does not improve for a given CPU time limit, stop
                    ctime=(cputime-tStart)/60;
                    if toDebug>=1
                        fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);fprintf(textfileName.interm_out,'Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);
                        disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                        fprintf('Termination criteria: Stop if fbest does not improve after %1.4f min. cputime holds.\n',IstopCondPara(11));
                    end
                    outputPara(7)=11;
                    break ; % exit while loop
               end
            end

            if IstopCondPara(4)~=0  % IstopCondPara(4)=1 check this condition, IstopCondPara(4)=0, skip it
                if (nIter - saveIter)> abs(IstopCondPara(4))  % If fbest does not improve after given num. of iter then stop the algo.
                    if (lastfbest - fbest)< 0.0001*lastfbest  % if the improvement in fbest is less than a tol, stop
                        ctime=(cputime-tStart)/60;
                        if toDebug>=1
                            fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);fprintf(textfileName.interm_out,'Final= %d; %f; %d; %g\n',nIter,fbest,nb,ctime);
                            disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                            fprintf('Termination criteria: Stop if fbest does not improve after %d iterations holds.\n',IstopCondPara(4));
                        end
                        outputPara(7)=4;
                        break ; % exit while loop
                    end
                end
            end

            if fbest<=(lastfbest-0.0001*lastfbest) % if fbest decreases by 0.01% of the last fbest then print the update
                saveIter=nIter;
                outputPara(11)=nIter;  % update the last fbest update iteration
                ctime=(cputime-tStart)/60;
                lastfbest=fbest;
                if toDebug>=1
                fprintf(textfileName.interm_out,'Interm.= %d; %f; %d; %g\n',nIter,fbest,nb,ctime); fprintf('Interm.= %d; %f; %d; %g\n',nIter,fbest,nb,ctime); 
                end
            end     

            
            % to save the results at an intermediate stage
            if toDebug>=1
                if IstopCondPara(4)<0 && intermFbestIterCnt<=nintermFbestIter
                    if isequal(nIter - saveIter,intermFbestIter(intermFbestIterCnt)), orgXstar=xbest.*normFactor; % get Xstar in the original space
                        fprintf(textfileName.interm_out,'intermFbestIter(%d)=%d, ',intermFbestIterCnt,(nIter - saveIter));
                        fprintf(textfileName.interm_out,'<%d> from position %d, ',nIter,deleteElt);fprintf(textfileName.interm_out,'CPU time=%f, ',ctime);fprintf(textfileName.interm_out,'numOfBoxesInList=%d, ',nb);fprintf(textfileName.interm_out,'fbest=%f, ',fbest);fprintf(textfileName.interm_out,'T.E.=%f, ',(norm(xMatrix*orgXstar-yArray))^2);fprintf(textfileName.interm_out,'P.E.=%f, ',((norm(xMatrix*orgXstar-xMatrix*trueb))^2)/((norm(xMatrix*trueb))^2));
                        fprintf(textfileName.interm_out,'Xstar=');printArray(orgXstar,'%d',textfileName.interm_out);fprintf(textfileName.interm_out,'\n');intermFbestIterCnt=intermFbestIterCnt+1;
                    end
                end
                if IstopCondPara(5)<0 && intermMaxIterCnt<=nintermMaxIter
                    if isequal( nIter,intermMaxIter(intermMaxIterCnt)),orgXstar=xbest.*normFactor; % get Xstar in the original space
                        fprintf(textfileName.interm_out,'intermMaxIter(%d)=%d, ',intermMaxIterCnt,nIter );
                        fprintf(textfileName.interm_out,'<%d> from position %d, ',nIter,deleteElt);fprintf(textfileName.interm_out,'CPU time=%f, ',ctime);fprintf(textfileName.interm_out,'numOfBoxesInList=%d, ',nb);fprintf(textfileName.interm_out,'fbest=%f, ',fbest);fprintf(textfileName.interm_out,'T.E.=%f, ',(norm(xMatrix*orgXstar-yArray))^2);fprintf(textfileName.interm_out,'P.E.=%f, ',((norm(xMatrix*orgXstar-xMatrix*trueb))^2)/((norm(xMatrix*trueb))^2));
                        fprintf(textfileName.interm_out,'Xstar=');printArray(orgXstar,'%d',textfileName.interm_out);fprintf(textfileName.interm_out,'\n');intermMaxIterCnt=intermMaxIterCnt+1;
                    end
                end
                if IstopCondPara(6)<0 && intermCPUCnt<=nintermCPU
                    if ctime>=intermCPU(intermCPUCnt), orgXstar=xbest.*normFactor; % get Xstar in the original space
                        fprintf(textfileName.interm_out,'intermCPU(%d)=%g, ',intermCPUCnt,ctime );
                        fprintf(textfileName.interm_out,'<%d> from position %d, ',nIter,deleteElt);fprintf(textfileName.interm_out,'CPU time=%f, ',ctime);fprintf(textfileName.interm_out,'numOfBoxesInList=%d, ',nb);fprintf(textfileName.interm_out,'fbest=%f, ',fbest);fprintf(textfileName.interm_out,'T.E.=%f, ',(norm(xMatrix*orgXstar-yArray))^2);fprintf(textfileName.interm_out,'P.E.=%f, ',((norm(xMatrix*orgXstar-xMatrix*trueb))^2)/((norm(xMatrix*trueb))^2));
                        fprintf(textfileName.interm_out,'Xstar=');printArray(orgXstar,'%d',textfileName.interm_out);fprintf(textfileName.interm_out,'\n');intermCPUCnt=intermCPUCnt+1;
                    end
                end    
            end

            % Select a new node to process and remove it from the list
            savecpu1=cputime;
            if secndSelOpt~=0  % if using secondary selection option as well
                if (nIter - saveIter)> switchAfterIter  % If fbest does not improve after a no. of iter. switch the selection criteria
                    if (nIter-preSwitchIter)> switchAfterIter 
                       if (lastfbest - fbest)< 0.0001*lastfbest % if the improvement in fbest is less than a tol, stop
                          if toDebug>=1
                             fprintf(textfileName.interm_out,'Selection criteria has been switched.\n');
                             fprintf('Selection criteria has been switched.\n');
                          end
                          switchSelOpt=1;
                       end
                       preSwitchIter=nIter;
                    end
                end
            end
            if switchSelOpt==1 % switch to secondary selection option
               switchSelOpt=0; 
               if secndSelOpt==4  % (BFS) box selection based on the lowest lb F(V)
                  selectnode = heap.extract_secondary(nBestBoxForSwitch, nRandBoxForSwitch, @(a,b) a.LBFBOX - b.LBFBOX );
                  isSelBFS=1;    % 1 if the selected box is using min lbF, =0 if any other selection option
               elseif secndSelOpt==5  % (DFS) select box with max #2 flags
                  selectnode = heap.extract_secondary(nBestBoxForSwitch, nRandBoxForSwitch, @(a,b) b.N2FLAG - a.N2FLAG ); 
                  isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
               elseif secndSelOpt==6  % select box with max #1 flags
                  selectnode = heap.extract_secondary(nBestBoxForSwitch, nRandBoxForSwitch, @(a,b) b.N1FLAG - a.N1FLAG );  
                  isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
               elseif secndSelOpt==7  % select box with max #0 flags
                  selectnode = heap.extract_secondary(nBestBoxForSwitch, nRandBoxForSwitch, @(a,b) b.N0FLAG - a.N0FLAG );
                  isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
               else
                  selectnode=heap.extract();  % buffer if any other secondary selection is defined 
                  isSelBFS=0;    % 1 if the selected box is using min lbF, =0 if any other selection option
               end  % secondary selection option
               
            else % use the primary selection option
                if firstSelOpt==1 % (BrFS) box selection based on age of the box criteria
                   selectnode=headnode;
                   if nb>1
                      headnode=headnode.Next; 
                   end
                   removeNode(selectnode);    % remove the selected node from the linked list
                else % other selection using Heap List
                   selectnode=heap.extract(); 
                end 
            end
            cpuArray(3)=cpuArray(3)+(cputime-savecpu1);
            nb=nb-1;

            % set selected node as a parent box
            Ybox=selectnode.BOX;               
            lbfY=selectnode.LBFBOX;             
            xlbfY=selectnode.LBFPT;             
            n2flagsY=selectnode.N2FLAG;         
            n0flagsY=selectnode.N0FLAG;     

            % special stopping condition if using BFS selection
            if fbest<=lbfY && isSelBFS==1  % if the selection is min lbF and fbest=lbf for the selected box, then stop
                if toDebug>=1
                   fprintf('Final= %d; %f; %d; %g\n',nIter,fbest,nb,(cputime-tStart)/60);disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$');
                   disp('Algo. 7 for BSS has stopped. fbest<=lbf for all the boxes in the list.');
                end
                outputPara(7)=0;
                break;  % exit while loop
            end
             
            % if using debugging mode
            if toDebug==2
               fprintf(textfileName.box_flag,'Iter=%d; deleteElt=%d; nb=%d \n',nIter,deleteElt,nb);
               fprintf(textfileName.box_flag,'\n');fprintf(textfileName.box_flag,'<%d> from position %d, ',nIter,deleteElt);
               fprintf(textfileName.box_flag,'parent box: ');printbox(n,Ybox,[],textfileName.box_flag);fprintf(textfileName.box_flag,'\n');
               fprintf(textfileName.box_flag,'lbf=%1.4f; ',lbfY);fprintf(textfileName.box_flag,'xlbf=');printArray(xlbfY,'%1.4f',textfileName.box_flag);fprintf(textfileName.box_flag,'\n');
            end
            fprintf(textfileName.lbFandfbest,'%1.8f; %1.8f \n',lbfY,fbest);
            
            savecpu1=cputime;
            %========================================================================================
            % bisect the selected box to get its child boxes and process them
            [lastfbestcpu,isStop_dfbestCPUHolds,necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,cpuArray,numOfBoxDel,nfEval,nFeval,nBoxDeleted,nb,xbest,fbest,lastnode,headnode,heap]=branch(heap,headnode,lastnode,n,nPts,tmax,Anorm,A_normdiv2,bnorm,cnorm,lbX,ubX,diagInv,xMatrix,yArray,Ybox,lbfY,xlbfY,n2flagsY,n0flagsY,decWidIdx,normxRelaxedOpt,cutDir,nCuts,cpuArray,...
                           textfileName,delCondPara,numOfBoxDel,nb,nBoxDeleted,nfEval,nFeval,toDebug,nIter,xbest,fbest,targetfbest,iPara,rPara,IotherPara,IstopCondPara,quadMinAlg,rankX,...
                           necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara,tStart,isStop_dfbestCPUHolds,lastfbestcpu);
            %====================================================================================================
            cpuArray(1)=cpuArray(1)+(cputime-savecpu1);    
            
        end  % end while stop==0

        
        totalcpu=cputime-tStart;
        outputPara(1)=nBoxDeleted;outputPara(2)=nb;outputPara(3)=cpuArray(3)/60;outputPara(4)=nFeval;outputPara(5)=nIter;outputPara(6)=cpuArray(6)/60;outputPara(9)=cpuIncFunEval/60;outputPara(12)=nfEval;
        outputPara(13)=cpuForRefineFesPt/60;outputPara(14)=nRefine;outputPara(15)=nFeval+nRefine;outputPara(16)=necConMaxVioRefQM;outputPara(17)=necConMaxVioInfQM;
        outputPara(18)=totalcpu;
        
       
         % save fbest ,xstar in a matlab data file to be used as warm start for the next run
         save(savedatafile,'fbest','xbest'); 
        
         if toDebug>=1
            fprintf(textfileName.interm_out,'IBB+, Set=%d, Total CPU =%1.8f min. for pDim=%d, tmax=%d \n',setCtr,totalcpu/60, n,tmax);
            fprintf(textfileName.interm_out,'CPU 0/0 for individual processess:======================== \n');
            fprintf(textfileName.interm_out,'(1) CPU : branch function =%1.8f  \n', (cpuArray(1))/(totalcpu) );
            fprintf(textfileName.interm_out,'(2) CPU : for expanding the list =%1.8f \n', cpuArray(2)/totalcpu);
            fprintf(textfileName.interm_out,'(3) CPU : selection and deletion of the box =%1.8f \n', (cpuArray(3))/(totalcpu) );
            fprintf(textfileName.interm_out,'(4) CPU : lb F eval =%1.8f \n', cpuIncFunEval/totalcpu);
            fprintf(textfileName.interm_out,'(5) CPU : fx eval =%1.8f \n', cpuArray(5)/totalcpu);
            fprintf(textfileName.interm_out,'(6) CPU : sampling =%1.8f \n', cpuArray(6)/totalcpu);
            fprintf(textfileName.interm_out,'(7) CPU : refinement =%1.8f \n', cpuForRefineFesPt/totalcpu);
            fprintf(textfileName.interm_out,'(8) CPU : bisecting the parent box =%1.8f \n', (cpuArray(8))/(totalcpu) );
            fprintf(textfileName.interm_out,'(9) CPU : process leaf nodes using DC2 and DC5 =%1.8f  \n', (cpuArray(9))/(totalcpu));
            fprintf(textfileName.interm_out,'(10) CPU : for selecting supp of feasible pt =%1.8f  \n', (cpuArray(4))/(totalcpu));
            fprintf(textfileName.interm_out,'(11) CPU : evaluate idx2flag in sampling =%1.8f  \n', (cpuArray(7))/(totalcpu));
            fprintf(textfileName.interm_out,'(12) cputime when last time fbest got updated =%1.8f min.  \n', (cpuArray(10))/(60));
            fprintf(textfileName.interm_out,'(2-9) sum of the above inidividual cpu =%1.8f min. \n', (cpuArray(1)+cpuArray(2)+cpuArray(3)+...
                cpuIncFunEval+cpuArray(5)+ cpuArray(6)+cpuForRefineFesPt+cpuArray(8)+cpuArray(9) )/60 );    
         end
         
        fprintf('IBB+, Set=%d, Total CPU =%1.8f min. for pDim=%d, tmax=%d \n',setCtr,totalcpu/60, n,tmax);
        fprintf('CPU 0/0 for individual processess:======================== \n');
        fprintf('(1) CPU : branch function =%1.8f  \n', (cpuArray(1))/(totalcpu) );
        fprintf('(2) CPU : for expanding the list =%1.8f \n', cpuArray(2)/totalcpu);
        fprintf('(3) CPU : selection and deletion of the box =%1.8f \n', (cpuArray(3))/(totalcpu) );
        fprintf('(4) CPU : lb F eval =%1.8f \n', cpuIncFunEval/totalcpu);
        fprintf('(5) CPU : fx eval =%1.8f \n', cpuArray(5)/totalcpu);
        fprintf('(6) CPU : sampling =%1.8f \n', cpuArray(6)/totalcpu);
        fprintf('(7) CPU : refinement =%1.8f \n', cpuForRefineFesPt/totalcpu);
        fprintf('(8) CPU : bisecting the parent box =%1.8f \n', (cpuArray(8))/(totalcpu) );
        fprintf('(9) CPU : process leaf nodes using DC2 and DC5 =%1.8f  \n', (cpuArray(9))/(totalcpu));
        fprintf('(10) CPU : for selecting supp of feasible pt =%1.8f  \n', (cpuArray(4))/(totalcpu));
        fprintf('(11) CPU : evaluate idx2flag in sampling =%1.8f  \n', (cpuArray(7))/(totalcpu));
        fprintf('(12) cputime when last time fbest got updated =%1.8f min.  \n', (cpuArray(10))/(60));
        fprintf('(2-9) sum of the above inidividual cpu =%1.8f min. \n', (cpuArray(1)+cpuArray(2)+cpuArray(3)+...
            cpuIncFunEval+cpuArray(5)+ cpuArray(6)+cpuForRefineFesPt+cpuArray(8)+cpuArray(9) )/60 ); 



%=======================================================================================================================================================
end  % end IBB
%=======================================================================================================================================================

%============================================================================================================================================================================    
function [lastfbestcpu,isStop_dfbestCPUHolds,necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,cpuArray,numOfBoxDel,nfEval,nFeval,nBoxDeleted,nb,xbest,fbest,lastnode,headnode,heap]=branch(heap,headnode,lastnode,n,nPts,tmax,Anorm,A_normdiv2,bnorm,cnorm,lbX,ubX,diagInv,xMatrix,yArray,Ybox,lbfY,xlbfY,n2flagsY,n0flagsY,decWidIdx,normxRelaxedOpt,cutDir,nCuts,cpuArray,...
                           textfileName,delCondPara,numOfBoxDel,nb,nBoxDeleted,nfEval,nFeval,toDebug,nIter,xbest,fbest,targetfbest,iPara,rPara,IotherPara,IstopCondPara,quadMinAlg,rankX,...
                           necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara,tStart,isStop_dfbestCPUHolds,lastfbestcpu)
     
        % get the child boxes after partitioning the parent box
        saveit=cputime;
        [Y_bisected,~,numOfChildBoxes,~,~]=bisectBoxAG72(n,Ybox,n2flagsY,n0flagsY,Anorm,bnorm,decWidIdx,normxRelaxedOpt,xlbfY,cutDir,nCuts);
        cpuArray(8)=cpuArray(8)+(cputime-saveit);

        % process all the child boxes --------------------------------------
        for j=1:numOfChildBoxes 

            isxbestupdated=0;  % flag to check if xbest gets updated using sampling
            if toDebug==2,fprintf(textfileName.box_flag,'child box: ');printbox(n,Y_bisected(j).box,[],textfileName.box_flag);end
            degIntvls=find(Y_bisected(j).box==0); 

            if IotherPara(3)~=0 
               if isStop_dfbestCPUHolds==1
                  return;  
               elseif ((cputime-lastfbestcpu)/60)>IstopCondPara(11)
                  isStop_dfbestCPUHolds=1;
                  return;  
               end
            end

            % check deletion conditions
            saveit=cputime;
            if delCondPara(2)==1  % DC2 
                if (Y_bisected(j).n2flags)>tmax % if #2 > k, the box is infeasible, delete the box
                   numOfBoxDel(2)=numOfBoxDel(2)+1;nBoxDeleted=nBoxDeleted+1;
                   if toDebug>=1
                      fprintf(textfileName.interm_out,'DC2.2-%d \n',nIter); if toDebug==2,fprintf(textfileName.box_flag,'DC2-%d fbest=%g \n',nIter,fbest);end
                   end
                   continue; % go the next iteration, which means discard the box

                elseif (Y_bisected(j).n2flags)==tmax   % if #2=k , use the point with lb Fbox as the feasible point, update fbest and discard the box                 
                   tempdegIntvls=find(Y_bisected(j).box~=2);
                   savecpu2=cputime; 
                   [necConMaxVioInfQM,VfesPt,ftp]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,xlbfY,diagInv,targetfbest,tempdegIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                   cpuIncFunEval=cpuIncFunEval+(cputime-savecpu2);
                   nFeval=nFeval+1;
                   if fbest > ftp
                      fbest=ftp;xbest=VfesPt;   %Ytilde=Y_bisected(j).box;  % save final solution
                      lastfbestcpu=cputime;cpuArray(10)=lastfbestcpu-tStart;
                   end
                   numOfBoxDel(2)=numOfBoxDel(2)+1;nBoxDeleted=nBoxDeleted+1;
                   if toDebug>=1
                      fprintf(textfileName.interm_out,'DC2.2-%d \n',nIter);if toDebug==2,fprintf(textfileName.box_flag,'DC2-%d fbest=%g \n',nIter,fbest);end
                   end
                   continue  % discard the box
                end
            end

            if delCondPara(5)==1   % DC5
                if ( n-Y_bisected(j).n0flags )==tmax   % if p-#0=k or #2+#1=k , use the point with lb Fbox as the feasible point, update fbest and discard the box
                    tempdegIntvls=find(Y_bisected(j).box==0); 
                    savecpu2=cputime;
                    [necConMaxVioInfQM,VfesPt,ftp]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,xlbfY,diagInv,targetfbest,tempdegIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                    cpuIncFunEval=cpuIncFunEval+(cputime-savecpu2);
                    nFeval=nFeval+1;
                    if fbest > ftp
                        fbest=ftp;xbest=VfesPt;   %Ytilde=Y_bisected(j).box;  % save final solution
                        lastfbestcpu=cputime;cpuArray(10)=lastfbestcpu-tStart;
                    end
                    numOfBoxDel(5)=numOfBoxDel(5)+1;nBoxDeleted=nBoxDeleted+1;
                    if toDebug>=1
                       fprintf(textfileName.interm_out,'DC5.0-%d \n',nIter);if toDebug==2,fprintf(textfileName.box_flag,'DC5-%d fbest=%g \n',nIter,fbest);end
                    end
                    continue  % discard the box
                end
            end
            
            if delCondPara(6)==1 % DC6, if #0 > p-k, then discard the box, redundant if no. of cuts is 1
                if IotherPara(22)==7 || IotherPara(22)==9 || IotherPara(22)==6 || IotherPara(22)==3 || IotherPara(22)==4 || IotherPara(22)==5
                    if Y_bisected(j).n0flags > (n-tmax)
                        numOfBoxDel(6)=numOfBoxDel(6)+1;nBoxDeleted=nBoxDeleted+1;
                        if toDebug>=1
                           fprintf(textfileName.interm_out,'DC6.0-%d \n',nIter);if toDebug==2,fprintf(textfileName.box_flag,'DC6-%d fbest=%g \n',nIter,fbest);end
                        end
                        continue;
                    end
                end
            end
            cpuArray(9)=cpuArray(9)+(cputime-saveit);

            % find a feasible point
            saveit=cputime;
            [necConMaxVioRefQM,nRefine,cpuForRefineFesPt,cpuArray,supp,VfesPt]=getFeasiblePt(n,degIntvls,find(Y_bisected(j).box~=2),lbX,ubX,xlbfY,targetfbest,tmax,IotherPara,IstopCondPara,Anorm,bnorm,cnorm,diagInv,...
                                                                    rPara,cpuArray,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara);
            cpuArray(6)=cpuArray(6)+(cputime-saveit);

            % update fbest if possible
            saveit=cputime;
            ftp=feval(VfesPt,A_normdiv2,bnorm,cnorm);
            cpuArray(5)=cpuArray(5)+(cputime-saveit);
            nfEval=nfEval+1;
            if fbest > ftp % if sampling found a better point, refine further using soft stop
                tempdegIntvls=1:n;tempdegIntvls(supp)=[];  %degIntvls=setdiff(1:n,supp);
                savecpu2=cputime;
                [necConMaxVioInfQM,VfesPt,ftp]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,VfesPt,diagInv,targetfbest,tempdegIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                cpuIncFunEval=cpuIncFunEval+(cputime-savecpu2);
                nFeval=nFeval+1;
                fbest=ftp;
                xbest=VfesPt;   %Ytilde=Y_bisected(j).box;  % save final solution
                isxbestupdated=1;
                lastfbestcpu=cputime;cpuArray(10)=lastfbestcpu-tStart;
            end

            % find lb f for the current box
            if Y_bisected(j).sdDP==1   % i.e. if the child box is the one with 0 flag, then find new lb f
                savecpu3=cputime;
                if IotherPara(7)==1
                    if (n-rankX)<Y_bisected(j).n0flags
                       [necConMaxVioInfQM,xlbf,v]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,xlbfY,diagInv,targetfbest,degIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                    else
                       v=lbfY;
                       xlbf=xlbfY; 
                    end
                else
                   [necConMaxVioInfQM,xlbf,v]=getLowerBound(n,Anorm,bnorm,cnorm,lbX,ubX,xlbfY,diagInv,targetfbest,degIntvls,iPara,rPara,quadMinAlg,1,necConMaxVioInfQM); % isSoftStop=1
                end
                cpuIncFunEval=cpuIncFunEval+(cputime-savecpu3);
                nFeval=nFeval+1;
            else   % the child box with '2' flag, we can use lb f from the parent box
                v=lbfY;
                xlbf=xlbfY;
            end

            % bound based deletion
            if delCondPara(1)==1  % if check DC1 is on
                if fbest <=v
                    nBoxDeleted=nBoxDeleted+1;
                    numOfBoxDel(1)=numOfBoxDel(1)+1;
                    if toDebug==2,fprintf(textfileName.box_flag,'DC1-%d lbFbox=%g fbest=%g \n',nIter,v,fbest);end
                    continue;
                end
            end % if delCondPara(1)==1

            % Call the branch function again, w/o adding the current box to the list
            if IotherPara(3)~=0 % then use this acceleration strategy
                if isxbestupdated==1
                    tempCutDir=zeros(1,tmax);
                    idxplus=1;
                    for ctr=1:tmax
                        if Y_bisected(j).box(supp(ctr))==1
                           tempCutDir(idxplus)=supp(ctr);
                           idxplus=idxplus+1;
                        end
                    end
                    tempCutDir=tempCutDir(1:idxplus-1); % remove the zeros
                    ntempCutDir=idxplus-1;
                    for ctr=1:min(ntempCutDir,IotherPara(3))
                       [lastfbestcpu,isStop_dfbestCPUHolds,necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,cpuArray,numOfBoxDel,nfEval,nFeval,nBoxDeleted,nb,xbest,fbest,lastnode,headnode,heap]=branch(heap,headnode,lastnode,n,nPts,tmax,Anorm,A_normdiv2,bnorm,cnorm,lbX,ubX,diagInv,xMatrix,yArray,Y_bisected(j).box,v,xlbf,Y_bisected(j).n2flags,Y_bisected(j).n0flags,decWidIdx,normxRelaxedOpt,tempCutDir(ctr),nCuts,cpuArray,...
                               textfileName,delCondPara,numOfBoxDel,nb,nBoxDeleted,nfEval,nFeval,toDebug,nIter,xbest,fbest,targetfbest,iPara,rPara,IotherPara,IstopCondPara,quadMinAlg,rankX,...
                               necConMaxVioInfQM,cpuIncFunEval,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara,tStart,isStop_dfbestCPUHolds,lastfbestcpu);
                    end
                end
            end

            % add the new boxes to the list
            
            nb=nb+1;
            saveit=cputime;
            if IotherPara(5)==1 % BrFS selection using Linked list
                newnode=dlnode(Y_bisected(j).box,v,xlbf,Y_bisected(j).n2flags,Y_bisected(j).n0flags,[],[]);   % create a new node 
                if nb==1  
                   headnode=newnode;
                else
                   newnode.insertAfter(lastnode);  % insert the new node after the last node
                end
                lastnode=newnode;                  % assign the newnode as last node now 
                clear newnode 
            else   
                heap.add( createbox(Y_bisected(j).box,v,xlbf,Y_bisected(j).n2flags,Y_bisected(j).n0flags,n-Y_bisected(j).n2flags-Y_bisected(j).n0flags,[],[]) );
            end
            cpuArray(2)=cpuArray(2)+(cputime-saveit); 
            if toDebug==2,fprintf(textfileName.box_flag,'lbFbox=%g fbest=%g \n',v,fbest); end

        end % for j=1:numOfChildBoxes

end  % function branch
      
%*********************************************************************************************************************

function  [necConMaxVioInfQM,xout,fout]=getLowerBound(n,Anorm,bnorm,cnorm,lb,ub,xtilde,diagInv,targetfbest,degIntvls,iPara,rPara,quadMinAlg,isSoftStop,necConMaxVioInfQM)

    isTF=0; % do not use targetfbest to find lb f inside quad. min. algo.
    
    newdim=n-length(degIntvls);  % dim of the reduced quadratic problem
    keepIdx=custom_setdiff(1:n,degIntvls); 
    Anew=Anorm(keepIdx,keepIdx); % initialize
    bnew=bnorm(keepIdx); % initialize
%         cnew=c_norm; % no change

    x0=xtilde( keepIdx );  % initial point for quad. min. algo.
    [Q,diagD]=eig(Anew,'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
    Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
    Sstruct.Q=Q(:, (newdim-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
    Sstruct.D=diagD( (newdim-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
    Sstruct.diagInv=diagInv(keepIdx);
    [epsMax,xquadmin,fout]=quad_min_algo_pool(newdim,Anew,bnew,cnorm,lb(keepIdx),ub(keepIdx),x0,Sstruct,iPara,rPara,quadMinAlg,isSoftStop,isTF,targetfbest);  % isSoftStop=1 means soft stop, isTF=0 do not use targetfbest  
    
    xout=zeros(n,1);xout(keepIdx)=xquadmin; % get solution x in the original space
    if necConMaxVioInfQM<epsMax, necConMaxVioInfQM=epsMax;end
  
    
end %***************************************************************************************************************


function [necConMaxVioRefQM,nRefine,cpuForRefineFesPt,cpuArray,supp,x_tilde]=getFeasiblePt(n,idx0flag,idx01flag,lbX,ubX,x0,targetfbest,tmax,IotherPara,IstopCondPara,...
                               A_norm,b_norm,c_norm,diagInv,rPara,cpuArray,cpuForRefineFesPt,nRefine,necConMaxVioRefQM,modify_iPara)
% x0 is a point such that f(x0)=lb f(V) for the working box V, or x0 will be xRelaxedOpt for the initial box X                                               
   isTF=0; % 0 means do not use the targetfbest inside quadratic minimization
   % Calculate the cardinality z of I={i:X_i contains 0}   
  
%    idx01flag=find(V~=2); % to save indices of the set i:X_i contains 0 
   saveit=cputime;
   idx2flag=custom_setdiff(1:n,idx01flag);
   cpuArray(7)=cpuArray(7)+(cputime-saveit);

   z=length(idx01flag); 
   
   tmax1=n-z;supp1=[];x_tilde=zeros(n,1);
   if tmax1>0
      [supp1,x_tilde]=feasiblePtL0_6(n,x0,tmax1,idx01flag,idx2flag,z,idx0flag);
   end
   tmax2=tmax-tmax1; idx0flag=custom_union(idx0flag,supp1);
   if tmax2>0
       idx12flag=1:n; idx12flag(idx0flag)=[];
       [~,localSupp,~,~,xreddim,~]=sfs(n-length(idx0flag),lbX(idx12flag),ubX(idx12flag),A_norm(idx12flag,idx12flag),b_norm(idx12flag),c_norm,tmax2,x0(idx12flag),[], diagInv(idx12flag),targetfbest,IotherPara,IstopCondPara,modify_iPara,rPara,[],0);  % toDebug=0,textfileName=[],ixstar=[]; 
       supp2=idx12flag(localSupp); % supp is global
       x_tilde(supp2)=xreddim(localSupp);
       supp=custom_union(supp1,supp2);
   else, supp=supp1;    
   end
  
   nRefine=nRefine+1;
   ctime=cputime; 
   nred=length(supp);   % reduce dim.
   [Q,diagD]=eig(A_norm(supp,supp),'vector');   % Q is an orthonormal matrix , diagD is a vector of eigenvalues in increasing order
   Sstruct.nNonZeroEig=sum(diagD>rPara(1));     % find no. of non zero eigenvalues
   Sstruct.Q=Q(:, (nred-Sstruct.nNonZeroEig+1):end );    % find Q matrix in the reduced space
   Sstruct.D=diagD( (nred-Sstruct.nNonZeroEig+1):end );  % find diagD in the reduced space
   Sstruct.diagInv=diagInv(supp);
   [epsMax,x_tilde(supp),~]=quad_min_algo_pool(nred,A_norm(supp,supp),b_norm(supp),c_norm,lbX(supp),ubX(supp),x_tilde(supp),Sstruct,modify_iPara,rPara,IotherPara(23),0,isTF,targetfbest);  % isSoftStop=0 hard stop, isTF=0 means do not use targetfbest 
   cpuForRefineFesPt=cpuForRefineFesPt+(cputime-ctime);
   if necConMaxVioRefQM<epsMax, necConMaxVioRefQM=epsMax;end
       
end

function [J,x_tilde]=feasiblePtL0_6(n,xrelax,k,I,J,z,degIntvls)
% 11 Feb22
% this subroutine will find the feasible point for L0 and perturbed L0 penalty for a given integer k
% n dim of the problem
% xrelax is the relaxedopt or OLS solution in the original dim
% z is the cardinality of I={i:X_i contains 0}
% k no. of non-zero comp we want to keep 
   
   % Step 3: 
   % 3a
%    J=setdiff(1:n,I);   % the set J={i:X_i does not contain 0}
%    J=custom_setdiff(1:n,I);   % the set J={i:X_i does not contain 0}
   
   %3b and 3c 
   x_tilde=zeros(n,1);

   x_tilde(J)=xrelax(J); % make the components in the set J non-zero
   small_vals=( abs(x_tilde(J))<eps ); % indices of small values
   x_tilde(small_vals)=sign(xrelax(small_vals)); % make those comp. 1 or -1 depending on the sign of the xrelax
        
   %3d     
   m=0;
   while m<(k-n+z)  % for k-n+z largest intervals
       
       [~,iv]=sort(abs(xrelax(I)),'descend'); % iv gives the indices of the largest intervals of I set in the decreasing order
       % find the index with the max width of X_i
       [id]=I(iv(1));

       if ~isempty(custom_intersect(degIntvls,id))   % if id is a member of degIntvls set
          I=custom_setdiff(I,id);
          continue ; 
       end
       
       if abs(xrelax(id))<eps   % if entry of xrelax < tol
          x_tilde(id)=sign(xrelax(id));  % keep it away from 0, 1 is an arbitrary no. any non-zero will work 
       else
          x_tilde(id)=xrelax(id); 
       end  
       I=custom_setdiff(I,id);
       J=custom_union(J,id);  % add id to the set J
       m=m+1;
   end

end

%**************************************************************************************************************

function [S,k,numOfChildBoxes,idx1flagY,nidx1flagY]=bisectBoxAG72(n,Ybox,n2flagsY,n0flagsY,A_norm,b_norm,decWidIdx,normxRelaxedOpt,xlbfY,cutDir,nCuts)
    
    idx1flagY=find( Ybox==1 );  
    nidx1flagY=n-n2flagsY-n0flagsY;
    if nidx1flagY<nCuts  % if the no. of cuts become more than the 1 flags boxes
        nCuts=nidx1flagY;
    end

    k=zeros(nCuts,1); % initialize, will store the components of cuts 
    if cutDir==-1  % choose the component i with max( width(X_i)) for i in I={j:X_j contains 0}
        
        save_width_1flags=intersect( decWidIdx,idx1flagY,'stable' );
%         save_width_1flags=custom_intersect( decWidIdx,idx1flagY );
        k=save_width_1flags(1:nCuts);  % take the first n_cuts entries as k

    elseif cutDir==-2 % choose the first 1 flag component
        k=idx1flagY(1:nCuts);


    elseif cutDir==-3  % choose the last 1 flag component

        k=idx1flagY( end-(nCuts-1):end);

    elseif cutDir==-4 % choose the comp which gives max univariate reduction in the quad term
        for j=1:nCuts 
%            k(j)=maxUnivarReduct1( setdiff(idx1flag,k(1:(j-1)),'stable') ,normxRelaxedOpt,A_norm);
           k(j)=maxUnivarReduct1( custom_setdiff(idx1flagY,k(1:(j-1)) ) ,normxRelaxedOpt,A_norm);
        end

    elseif cutDir==-5 
        for j=1:nCuts 
%            k(j)=maxUnivarReduct2( setdiff(idx1flag,k(1:(j-1)),'stable') ,A_norm,b_norm);
           k(j)=maxUnivarReduct2( custom_setdiff(idx1flagY,k(1:(j-1)) ) ,A_norm,b_norm);
        end

    elseif cutDir==-6
        for j=1:nCuts 
%            k(j)=minUnivarReduct1( setdiff(idx1flag,k(1:(j-1)),'stable') ,normxRelaxedOpt,A_norm);
           k(j)=minUnivarReduct1( custom_setdiff( idx1flagY,k(1:(j-1)) ) ,normxRelaxedOpt,A_norm);
        end
    
    elseif cutDir==-8
       [~,k]=maxk(abs(xlbfY(idx1flagY)),nCuts);
       k=idx1flagY(k);

    else % partition the direction given by cutDir
       k=cutDir; 

    end

    indvec=zeros(nCuts,1);
    S=struct('sdDP',[],'box',[],'n2flags',[],'n0flags',[]);
    [S,~,~]=getAllVertices(nCuts,nCuts,k,Ybox,n2flagsY,n0flagsY,nCuts,indvec,S);
    S(1)=[]; % delete the first empty box in the structure
    numOfChildBoxes=2^(nCuts);


end

%*****************************************************************************************************************************

function [S,n0flagsX,n2flagsX]=getAllVertices(ii,n,k,X,n2flagsX,n0flagsX,n_cuts,indvec,S)
% ii=length(k)
% n=length(k)
% k=indices of the components to be cut
% I interval vector of boxes with lower corner
% J interval vector of boxes with upper corner
% X original box 
% indvec initialization of nx1 scalar vector
% Z is empty structure to save the output boxes

    for i=1:2
        indvec(ii)=i;
        if ii>1
        [S,n0flagsX,n2flagsX]=getAllVertices(ii-1,n,k,X,n2flagsX,n0flagsX,n_cuts,indvec,S);
        else
          
            boxflagtemp=X;
            n2flags=n2flagsX;n0flags=n0flagsX;
            for j=1:n_cuts
               if indvec(j)==1
                  sdDPflag=1;
                  boxflagtemp(k(j))=0;
                  n0flags=n0flags+1;
                     
               elseif indvec(j)==2
                  if sum(indvec)==(2*n_cuts)
                     sdDPflag=0; 
                  else
                     sdDPflag=1;
                  end
                  boxflagtemp(k(j))=2;
                  n2flags=n2flags+1;
                 
               end
            end
            len=length(S)+1;
            S(len).box=boxflagtemp;
            S(len).sdDP=sdDPflag;
            S(len).n2flags=n2flags;
            S(len).n0flags=n0flags;
 
        end
        
    end
end

%*************************************************************************************************************

function [d]=maxUnivarReduct1(supp,normxRelaxedOpt,A)
  
  %1. direction to partition the box will be among those which are in the
  %support of the box
  nc=length(supp); % number of components
  A=A(supp,supp);  % A in the reduced dim
  xRelaxOpt=normxRelaxedOpt(supp); % only those comp which are in the support  (possible improvement, get the ols sol again in the reduced dim "maxUnivarReduct2")
  xref=xRelaxOpt'*A*xRelaxOpt;  % the reference value to which we will compare the reduction
  reduct=zeros(1,nc); % to save the reduction due to each component
  for i=1:nc
      xp_minus_1=xRelaxOpt;xp_minus_1(i)=0;  % after removing one comp
      reduct(i)=xref - (xp_minus_1'*A*xp_minus_1);  % save each univariate reduction      
  end
  [~,idx]=max(reduct); % the comp that gives the most reduction is chosen
  d=supp(idx);
  
  
end

%***********************************************************************************************************

function [d]=maxUnivarReduct2(supp,A,b)
% Difference bw maxUnivarReduct1 and this is that here we are again computing xols in the reduced dim determined by the support


  %1. direction to partition the box will be among those which are in the
  %support of the box
  nc=length(supp); % number of components
  A=A(supp,supp);  % A in the reduced dim
  b=b(supp);  % b in the reduced dim
  xRelaxOpt=A\(-b); % only those comp which are in the support  (possible improvement, get the ols sol again in the reduced dim ???)
  xref=xRelaxOpt'*A*xRelaxOpt;  % the reference value to which we will compare the reduction
  reduct=zeros(1,nc); % to save the reduction due to each component
  for i=1:nc
      xp_minus_1=xRelaxOpt;xp_minus_1(i)=0;  % after removing one comp
      reduct(i)=xref - (xp_minus_1'*A*xp_minus_1);  % save each univariate reduction      
  end
  [~,idx]=max(reduct); % the comp that gives the most reduction is chosen
  d=supp(idx);
  
  
end

%**************************************************************************************************************

function [d]=minUnivarReduct1(supp,normxRelaxedOpt,A)
  
  %1. direction to partition the box will be among those which are in the
  %support of the box
  nc=length(supp); % number of components
  A=A(supp,supp);  % A in the reduced dim
  xRelaxOpt=normxRelaxedOpt(supp); % only those comp which are in the support  (possible improvement, get the ols sol again in the reduced dim "maxUnivarReduct2")
  xref=xRelaxOpt'*A*xRelaxOpt;  % the reference value to which we will compare the reduction
  reduct=zeros(1,nc); % to save the reduction due to each component
  for i=1:nc
      xp_minus_1=xRelaxOpt;xp_minus_1(i)=0;  % after removing one comp
      reduct(i)=xref - (xp_minus_1'*A*xp_minus_1);  % save each univariate reduction      
  end
  [~,idx]=min(reduct); % the comp that gives the most reduction is chosen
  d=supp(idx);
  
  
end %********************************************************************************************************************

function printbox(n,X,numOfBoxCtr,textfileid) 
   
   if isempty(numOfBoxCtr)
      string='[%g'; 
   else
      string='%g[%g';
   end
   for i=2:(n-1)
      string=[string '%g'];
   end
   string=[string '%g]  '];
   fprintf(textfileid,string,[numOfBoxCtr X]);
   
    
end %****************************************************************************************************************************

function [out]=feval(x,Adiv2,b,c)  
%{
Input= x is a floating number
Output= a floating number f(x)
yArray    % n-dim row vector of given y values
xMatrix   % pxn matrix from the given b_i points
%}


  % if LLS or Unconst. Algo.
  out=x'*(Adiv2*x+b);
  out=out+c;
      
end

%************************************************************************************************************************************

function Z = custom_setdiff(X,Y)
% https://www.mathworks.com/matlabcentral/answers/53796-speed-up-intersect-setdiff-functions#comment_1002187
% Idea from-  Afaik Kevin Murphy of the Bayes Net TB

    if ~isempty(X)&&~isempty(Y)
      check = false(1, max(max(X), max(Y)));
      check(X) = true;
      check(Y) = false;
      Z = X(check(X));  
    else
      Z = X;
    end
end

%************************************************************************************************************************************

function C = custom_intersect(A,B)
% https://www.mathworks.com/matlabcentral/answers/53796-speed-up-intersect-setdiff-functions#comment_1002187
% Idea from-  Afaik Kevin Murphy of the Bayes Net TB

    if ~isempty(A)&&~isempty(B)
     P = zeros(1, max(max(A),max(B)) ) ;
     P(A) = 1;
     C = B(logical(P(B)));
    else
      C = [];
    end

end

%************************************************************************************************************************************

function U=custom_union(A,B)
  
   if isempty(A)
      U=B;
   elseif isempty(B)
      U=A;
   else
      check = zeros(1, max(max(A), max(B)));
      check(A) = A;
      check(B) = B;
      U = check(logical(check)); 
   end

end %*****************************************************************************************************************************
