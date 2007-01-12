       subroutine permutation(nperm,nrow,ncol,ngeno,model,
     +                        genotRate,data,Gstat)
       
       implicit none
       integer nperm,nrow,ncol,ngeno,genotRate,model
       integer data(nrow,ncol)
       double precision Gstat(ncol-1,nperm)
       integer i,j
       integer p(nrow),aux(nrow)
       
       do i=1,nrow
        aux(i)=data(i,1) 
       end do    
 
c initiate the n permutations (nperm)
       do i=1,nperm
c permutation of case-control label
        call rperm(p,nrow) 
        do j=1,nrow
         data(j,1)=aux(p(j)) 
        end do
c computing G statistic
        call WGassociation(nrow,ncol,ngeno,model,genotRate,
     +          data,Gstat(:,i))
       end do

       end subroutine





       subroutine WGassociation(nrow,ncol,ngeno,model,genotRate,
     +                         data,Gstat)

       implicit none
       integer nrow,ncol,ngeno,genotRate,ipred,i,j
       integer data(nrow,ncol),nn(2,ngeno),model
       double precision Gstat(ncol-1)
       real x(nrow), y(nrow),control, r2
       integer nomis
    
       do i=2,ncol
        if (model.ne.5) then
         call table(nrow,ngeno,data(:,1),data(:,i),model,genotRate,
     +            nn,ipred)
         if (ipred.eq.1) then
          call G(ngeno,nn,Gstat(i-1))
         else
          Gstat(i-1)=-2         
         end if

c additive
        else if (model.eq.5) then
         nomis=0
         do j=1,nrow
           if (data(j,i).ne.0) then
            x(nomis+1)=real(data(j,1))
            y(nomis+1)=real(data(j,i))
            nomis=nomis+1
           end if  
          end do

         control=(real(nomis)/real(nrow))*100.0
         if (control.lt.genotRate) then
          Gstat(i-1)=-2
         else if (minval(y(1:nomis)).eq.maxval(y(1:nomis))) then
          Gstat(i-1)=-1
         else
          
c VM
          Gstat(i-1)=(nomis-1)*r2(x,y,nomis)
         end if
        end if  

       end do 
      
       end subroutine

    



      subroutine table(n,geno,vec1,vec2,model,genotRate,nn,ipred)
      implicit none

c
c   IN:
c   n -- length of vectors (number of individuals)
c   geno -- 2 or 3 indicating the number of genotypes (to collapse dominant, recessive, etc.)
c   vec1 -- 1 y 0
c   vec2 -- 1,2 or 1,2,3
c   genotRate -- percentage of desired genotyping 
c
c   OUT:
c   nn -- 2x3 1:case  2:control  
c   ipred -- controls wheter the genotyping is less than genotRate
c
   
      integer n,i,j,k,geno,genotRate,ipred,nOK
      integer vec1(n),vec2(n)
      integer nn(2,geno),model
      real control

      do j=1,2
       do k=1,geno
        nn(j,k)=0
       end do
      end do          


c  Codominant
      if (model.eq.1) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.1) then
            nn(1,1)=nn(1,1)+1
         else if (vec2(i).eq.2) then
            nn(1,2)=nn(1,2)+1
         else if (vec2(i).eq.3) then
            nn(1,3)=nn(1,3)+1 
         end if
        else
         if (vec2(i).eq.1) then
            nn(2,1)=nn(2,1)+1
         else if (vec2(i).eq.2) then
            nn(2,2)=nn(2,2)+1
         else if (vec2(i).eq.3) then
            nn(2,3)=nn(2,3)+1 
         end if
        end if        
       end do   

c  Dominant
      else if (model.eq.2) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.1) then
            nn(1,1)=nn(1,1)+1
         else if ((vec2(i).eq.2).or.(vec2(i).eq.3)) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if (vec2(i).eq.1) then
            nn(2,1)=nn(2,1)+1
         else if ((vec2(i).eq.2).or.(vec2(i).eq.3)) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

c  Recessive
      else if (model.eq.3) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if ((vec2(i).eq.1).or.(vec2(i).eq.2)) then
            nn(1,1)=nn(1,1)+1
         else if (vec2(i).eq.3) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if ((vec2(i).eq.1).or.(vec2(i).eq.2)) then
            nn(2,1)=nn(2,1)+1
         else if (vec2(i).eq.3) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

c  Overdominant
      else if (model.eq.4) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.2) then
            nn(1,1)=nn(1,1)+1
         else if ((vec2(i).eq.1).or.(vec2(i).eq.3)) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if (vec2(i).eq.2) then
            nn(2,1)=nn(2,1)+1
         else if ((vec2(i).eq.1).or.(vec2(i).eq.3)) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

      end if



      nOK=0
      do j=1,2
       do k=1,geno
        nOK=nOK+nn(j,k)
       end do
      end do          
   
      ipred=1
      control=(real(nOK)/real(n))*100.0
      if (control.lt.genotRate) then
        ipred=0
      end if
  
      end subroutine table 


      subroutine G(geno,tt,Gstat)
      implicit none
      integer geno
      integer tt(2,geno),i,j
      double precision r(2),s(geno),N,Gstat,Gij,obs,esp

      do i=1,2
       if (geno.eq.3) then 
        r(i)=tt(i,1)+tt(i,2)+tt(i,3)
       else if (geno.eq.2) then
        r(i)=tt(i,1)+tt(i,2)
       end if 
      end do
      N=r(1)+r(2)

      do i=1,geno
       s(i)=tt(1,i)+tt(2,i)
      end do


     
      if ((geno.eq.3).and.(((s(1).eq.0).and.(s(2).eq.0)).or.
     +   ((s(1).eq.0).and.(s(3).eq.0)).or.
     +   ((s(2).eq.0).and.(s(3).eq.0)))) then

         Gstat=-1.0

      else if ((geno.eq.2).and.((s(1).eq.0).or.(s(2).eq.0).or.
     +   (s(3).eq.0))) then
    
         Gstat=-1.0

      else
       Gstat=0.0
       do i=1,2
        do j=1,geno
         if (tt(i,j).gt.0) then
          esp=(r(i)*s(j))/N
          obs=tt(i,j)
          Gij=obs*log(obs/esp)
          Gstat=Gstat+Gij  
         end if
         end do
       end do   
       Gstat=2.0*Gstat  
      end if  

      end subroutine G 




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      SUBROUTINE RPERM(P,N)
C===Generate a random permutation, P, of the first N integers.
C   (equivalent to sampling WITHOUT REPLACEMENT).
C   Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
      INTEGER N,P(N), K,J,I,IPJ,ITEMP,M
      REAL U(100)
      DO 1 I=1,N
1     P(I)=I
C---Generate up to 100 U(0,1) numbers at a time.
      DO 3 I=1,N,100
        M=MIN(N-I+1,100)
        CALL RAND(U,M)
        DO 2 J=1,M
          IPJ=I+J-1
          K=INT(U(J)*(N-IPJ+1))+IPJ
          ITEMP=P(IPJ)
          P(IPJ)=P(K)
2       P(K)=ITEMP
3     CONTINUE
      RETURN
      END
      SUBROUTINE IRAND(S,N,LOW,HI)
C===Generate a random integer sequence: S(1),S(2), ... ,S(N)
C   such that each element is in the closed interval <LOW,HI> and
C   sampled WITH REPLACEMENT.                            HDK, JUNE 1971.
      INTEGER N,S(N),LOW,HI,IX,I
      REAL U(1)
      DOUBLE PRECISION X
      DO 1 I=1,N
        CALL RAND(U,1)
C---Use DP arithmetic to effect a more precise transformation.
        X=DBLE((HI+1)-LOW)*U(1) + DBLE(LOW)
        IX=X
        IF(X.LT.0 .AND. IX.NE.X) IX=X-1.D0
        S(I)=IX
1     CONTINUE
      RETURN
      END
      BLOCK DATA
C=======================================================================
C  Portable pseudo-random integer generator, especially for
C  microcomputers with a limitation of 16 bit integers. Translated from
C  original Pascal version(1) to Fortran 77 by H. D. Knoble, PSU.
C
C   The supporting paper is:
C   (1) B. Wichmann & D. Hill, "Building a Random-Number Generator",
C             BYTE, March, 1987, 12(3):127-128.
C
C   Also see the following related works:
C   (2) Wichmann, B.A. and Hill, I.D., "An Efficient and Portable",
C             Pseudo-random Number Generator", Applied Statistics,
C             Volume 31, 1982, pages 188-190.
C   (3) Haas, Alexander, "The Multiple Prime Random Number Generator",
C             ACM Transactions on Mathematical Software; December,
C             1987, 13(4):368-381.
C   (4) L'Ecuyer, Pierre, "Efficient and Portable Combined Random Number
C             Generators", Communications of the ACM; June, 1988,
C             31(6):742-749,774.
C
C Use...
C      CALL RAND(U,N)
C          To generate a sequence, U, of N Uniform(0,1) numbers.
C          Cycle length is ((30269-1)*(30307-1)*(30323-1))/4  or
C          6953607871644  > 6.95E+12.
C
C     To access the SEED vector in the calling program use statements:
C     INTEGER SEED(3)
C     COMMON/RANDOM/SEED
C
C  The common variable SEED is the array of three current seeds.
      INTEGER SEED(3)
      COMMON/RANDOM/SEED
      DATA SEED(1),SEED(2),SEED(3)/1,10000,3000/
      END
C=======================================================================
      SUBROUTINE RAND(U,N)
      INTEGER N,X,Y,Z
      REAL U(N),V
      COMMON/RANDOM/X,Y,Z
      IF(N.LE.0) RETURN
      DO 1 I=1,N
        X=171*MOD(X,177)-2*(X/177)
        IF(X.LT.0) X=X+30269
        Y=172*MOD(Y,176)-35*(Y/176)
        IF(Y.LT.0) Y=Y+30307
        Z=170*MOD(Z,178)-63*(Z/178)
        IF(Z.LT.0) Z=Z+30323
        V=X/30269.0 + Y/30308.0 + Z/30323.0
1       U(I)=V-INT(V)
      RETURN
      END



      REAL FUNCTION r2(x,y,n)

      INTEGER n
      REAL x(n),y(n)

      INTEGER i
      REAL sxx,sxy,syy,sx,sy, nn
      nn=real(n)

      sx=0.
      sy=0.
      sxx=0.
      syy=0.
      sxy=0.
      do 12 i=1,n
        sx=sx+x(i)
        sy=sy+y(i)
        sxx=sxx+x(i)*x(i)
        syy=syy+y(i)*y(i)
        sxy=sxy+x(i)*y(i)
12    continue

      r2=((sxy-sx*sy/nn)**2)/(sxx-sx*sx/nn)/(syy-sy*sy/nn)

      return
      END


