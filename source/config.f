      module config
      integer nmax,msrs,nmaxl,mmax,nspscale,
     1        nspirals,emax,maxit
      PARAMETER (Nmax=5000)    !  Change NMAX if RAM is not sufficient
      PARAMETER (msrs=56+1)     !  Size of Schlegel output matrix
      PARAMETER (NmaxL=(Nmax*(Nmax-1))/2)
      PARAMETER (Mmax=Nmax/2+2)
      PARAMETER (NSpScale=12)
      PARAMETER (NSpirals=Mmax*NSpScale)
      PARAMETER (Emax=3*Nmax/2)
      PARAMETER (maxit=2000000)       
      save
      end module config
