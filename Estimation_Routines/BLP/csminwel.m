function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
% [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
% Clean/silent version: suppresses all output unless Verbose=1.

[nx,no]=size(x0);
nx=max(nx,no);
Verbose=0; % Set to 1 for output, 0 for silent
NumGrad= isempty(grad);
done=0;
itct=0;
fcount=0;

f0 = feval(fcn,x0,varargin{:});
if f0 > 1e50
    if Verbose, disp('Bad initial parameter.'), end
    return
end
if NumGrad
   if length(grad)==0
      [g badg] = numgrad(fcn,x0, varargin{:});
   else
      badg=any(find(grad==0));
      g=grad;
   end
else
   [g badg] = feval(grad,x0,varargin{:});
end

retcode3=101;
x=x0;
f=f0;
H=H0;

while ~done
   g1=[]; g2=[]; g3=[];
   itct=itct+1;
   [f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,varargin{:});
   fcount = fcount+fc;

   if retcode1 ~= 1
      if retcode1==2 || retcode1==4
         wall1=1; badg1=1;
      else
         if NumGrad
            [g1 badg1] = numgrad(fcn, x1,varargin{:});
         else
            [g1 badg1] = feval(grad,x1,varargin{:});
         end
         wall1=badg1;
         % Debugging save; rarely needed
         if Verbose
             % save('g1.mat', 'g1', 'x1', 'f1', 'varargin');
         end
      end
      if wall1 && (length(H) > 1)
         Hcliff=H+diag(diag(H).*rand(nx,1));
         if Verbose, disp('Cliff.  Perturbing search direction.'), end
         [f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,varargin{:});
         fcount = fcount+fc;
         if  f2 < f
            if retcode2==2 || retcode2==4
                  wall2=1; badg2=1;
            else
               if NumGrad
                  [g2 badg2] = numgrad(fcn, x2,varargin{:});
               else
                  [g2 badg2] = feval(grad,x2,varargin{:});
               end
               wall2=badg2;
               if Verbose
                  % badg2
                  % save g2 g2 x2 f2 varargin
               end
            end
            if wall2
               if Verbose, disp('Cliff again.  Try traversing'), end
               if norm(x2-x1) < 1e-13
                  f3=f; x3=x; badg3=1;retcode3=101;
               else
                  gcliff=((f2-f1)/((norm(x2-x1))^2))*(x2-x1);
                  if(size(x0,2)>1), gcliff=gcliff'; end
                  [f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),varargin{:});
                  fcount = fcount+fc;
                  if retcode3==2 || retcode3==4
                     wall3=1; badg3=1;
                  else
                     if NumGrad
                        [g3 badg3] = numgrad(fcn, x3,varargin{:});
                     else
                        [g3 badg3] = feval(grad,x3,varargin{:});
                     end
                     wall3=badg3;
                     if Verbose
                         % badg3
                         % save g3 g3 x3 f3 varargin;
                     end
                  end
               end
            else
               f3=f; x3=x; badg3=1; retcode3=101;
            end
         else
            f3=f; x3=x; badg3=1;retcode3=101;
         end
      else
         f2=f; f3=f; badg2=1; badg3=1; retcode2=101; retcode3=101;
      end
   else 
      f2=f;f3=f;f1=f;retcode2=retcode1;retcode3=retcode1;
   end
   %how to pick gh and xh
   if f3 < f - crit && badg3==0
      %ih=3
      fh=f3;xh=x3;gh=g3;badgh=badg3;retcodeh=retcode3;
   elseif f2 < f - crit && badg2==0
      %ih=2
      fh=f2;xh=x2;gh=g2;badgh=badg2;retcodeh=retcode2;
   elseif f1 < f - crit && badg1==0
      %ih=1
      fh=f1;xh=x1;gh=g1;badgh=badg1;retcodeh=retcode1;
   else
      [fh,ih] = min([f1,f2,f3]);
      if Verbose, disp(['ih = ',num2str(ih)]); end
      switch ih
         case 1
            xh=x1;
         case 2
            xh=x2;
         case 3
            xh=x3;
      end
      retcodei=[retcode1,retcode2,retcode3];
      retcodeh=retcodei(ih);
      if exist('gh','var')
         nogh=isempty(gh);
      else
         nogh=1;
      end
      if nogh
         if NumGrad
            [gh badgh] = numgrad(fcn,xh,varargin{:});
         else
            [gh badgh] = feval(grad, xh,varargin{:});
         end
      end
      badgh=1;
   end
   stuck = (abs(fh-f) < crit);
   if (~badg)&&(~badgh)&&(~stuck)
      H = bfgsi(H,gh-g,xh-x);
   end
   if Verbose
      disp('----');
      disp(sprintf('Improvement on iteration %d = %18.9f',itct,f-fh));
      if itct > nit
         disp('iteration count termination');
      elseif stuck
         disp('improvement < crit termination');
      end
      rc=retcodeh;
      if rc == 1
         disp('zero gradient');
      elseif rc == 6
         disp('smallest step still improving too slow, reversed gradient');
      elseif rc == 5
         disp('largest step still improving too fast');
      elseif (rc == 4) || (rc==2)
         disp('back and forth on step length never finished');
      elseif rc == 3
         disp('smallest step still improving too slow');
      elseif rc == 7
         disp('warning: possible inaccuracy in H matrix');
      end
   end
   f=fh;
   x=xh;
   g=gh;
   badg=badgh;
   if itct > nit || stuck
       done = 1;
   end
end
end
