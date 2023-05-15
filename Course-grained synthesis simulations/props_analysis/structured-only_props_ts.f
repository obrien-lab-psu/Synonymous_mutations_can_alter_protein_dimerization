*--------------------------------------------------------------*
      PROGRAM  MODTRAJ
*
*  This program modifies the binary charmm trajectory header
*  Author: ANONYMOUS
*
*  Read in file.
*
*  Contact types:
*    Type I contacts are between residues within the same 2' element
*    Type II contacts are between residues in different 2' elements
*    Type III contacts are between one residue in a 2' element and one
*    residue that is not in a 2' element
*    Type IV contacts are between two residues not in 2' elements;
*    these contacts are not counted in the current version
*
* Declare some variables
      IMPLICIT NONE
      CHARACTER*128 INFILE, OUTFILE, NSREF, CRD, model, CNTRL
      CHARACTER*4 hdrr
      CHARACTER*80 line(2),linecrd(10)
      INTEGER nat,nres,ntot,sc(500),scat(500),cnres
      INTEGER nskip,bb(500),bbat(500),nfreat,nsecstruc
      INTEGER cnatBB,cnatBS,cnatSS,contactmapBB(500,500)
      INTEGER contactmapBS(500,500),contactmapSS(500,500),nstruc
      INTEGER cinatBB,cinatBS,cinatSS,cinatsbb(1000),cnatsbb(1000)
      INTEGER cinnatBB,cinnatBS,cinnatSS,cinattot,cnattot,cinnattot
      INTEGER typeonetot, typetwotot, typethreetot, use_num
      INTEGER typeonei, typetwoi, typethreei
      INTEGER I,J,K,L,m,natrec,nfret,nframes,struc(100,100),val(100),f
      INTEGER*4 icntrl4(20),natrec4
      INTEGER*4 ntitle
      INTEGER*8 jcntrl8(20),icntrl8(20)
      REAL*4 X(500),Y(500),Z(500),Rcm(3),temp,dist,natcut,Rg,RgNS,e2ed
      REAL*4 distmapBS(500,500),distmapSS(500,500),pij(500,500)
      REAL*4 distmapBB(500,500),maxdist,sdist
      CHARACTER(len=60) :: arg

* get arguements from command line
      DO i = 1, iargc()
        CALL getarg(i, arg)
        if(i.eq.1) NSREF=arg
        if(i.eq.2) CRD=arg
        if(i.eq.3) CNTRL=arg
        if(i.eq.4) INFILE=arg
        if(i.eq.5) OUTFILE=arg
        if(i.eq.6) then
* convert string to real
         read (arg,*) temp
* and real to integer
         nsecstruc=int(temp)
        endif
        if(i.eq.7) then
* convert string to real
         read (arg,*) temp
* and real to integer
         nres=int(temp)
        endif
        if(i.eq.8) then
* convert string to real
         read (arg,*) temp
* and real to integer
         ntot=int(temp)
        endif
        if(i.eq.9) then
* arg must either be 'ca' or 'cacb'         
         model=arg
        endif
        if(i.gt.9) then
         print *, 'ERROR: only allowed 8 commands'
         STOP
        endif
        print *, i, arg
      END DO
* define some hardcoded parameters
* number of atoms
      nat = ntot
* number of lines to skip in CRD file
      nskip = 6
* for contact map
      NATCUT = 8
* scale factor for defining distance cutoff of native contacts
      sdist = 1.2 
      if(NSREF .ne. 'none') then
*
* open native state structure file, in xyz format
      OPEN(UNIT=55,FILE=NSREF,
     $      FORM='FORMATTED')
      endif
*
* open input trajectory file 
      OPEN(UNIT=56,FILE=INFILE,FORM='UNFORMATTED')

* open output data file
      OPEN(57,FILE=OUTFILE, STATUS='UNKNOWN',
     $      FORM='formatted')

* open CRD file to ID backbone and sidechain atoms
      OPEN(UNIT=58,FILE=CRD,
     $      FORM='FORMATTED')

* set all residues to default of no sidechain atoms
      DO I=1,nres
        sc(I)=0
        scat(I)=0
        PRINT*,I,nres
      ENDDO
      cnres = 0
      nskip = 2

      DO I=1,nat+nskip
        IF(I.LT.nskip+1) THEN
          READ(58,*)
        ELSE
          READ(58,*) (linecrd(J),J=1,9)
          IF((linecrd(4) .eq. 'B').and.(model .eq. 'cacb')) THEN
            sc(cnres) = 1
            scat(cnres) = I-nskip
          ELSEIF ((linecrd(4) .eq. 'A').and.(model .eq. 'cacb')) THEN
            cnres = cnres+1
            bb(cnres) = 1
            bbat(cnres) = I-nskip
          ELSEIF ((linecrd(4) .eq. 'A').and.(model .eq. 'ca')) THEN
            cnres = cnres+1
            bb(cnres) = 1
            bbat(cnres) = I-nskip
          ELSE
            PRINT*, 'ERROR: Unrecognized atom type', linecrd(4)
            STOP
          ENDIF
        ENDIF
      ENDDO
      IF(cnres .NE. nres) THEN
        print*, cnres,'NE',nres
        print*, 'ERROR: incorrect number of residues'
        STOP
      ENDIF
      DO I=1,nres
        PRINT*, I, bbat(I), bb(I)
        if(model .eq. 'cacb') PRINT*, I, scat(I), sc(I)
      ENDDO
      close(58)

* read in residues for each secondary structure
      print*,'Secondary Structure residues'
      OPEN(UNIT=58,FILE=CNTRL,
     $      FORM='FORMATTED')
       do j=1,nsecstruc
        READ(58,*) (val(i),i=1,3)
        struc(val(1),1)=val(2)
        struc(val(1),2)=val(3)
        nstruc=val(1)
        PRINT*, nstruc, struc(val(1),1), struc(val(1),2)
       enddo
      close(58)

* read in header and control info as i*4
      READ (56) hdrr,icntrl4
* convert the i*4 to i*8
      DO I=1,20
        jcntrl8(I) = icntrl4(I)
        icntrl8(I) = icntrl4(I)
      ENDDO
      PRINT *, hdrr, icntrl4
* just set the first step of every file to 0
      icntrl4(2) = 0

* read in the number of title lines
      READ (56) ntitle, (line(I),I=1,ntitle)
      PRINT *, ntitle
      DO I=1,ntitle
        PRINT *, 'HERE', line(I)
      ENDDO
* read in the number of title line
      READ (56) natrec4
      PRINT *, natrec4
      natrec = natrec4
* do a quality control check
      NFREAT=NATREC-ICNTRL8(9)
      IF(NFREAT.NE.NATREC) THEN
        PRINT*, 'ERROR: NFREAT .NE. NATREC'
        !STOP
      ENDIF

      if(NSREF .ne. 'none') then
* read in coordinates of native state,
* and compute its distance matrix
       READ(55,*) (X(i),i=1,ntot)
       READ(55,*) (Y(i),i=1,ntot)
       READ(55,*) (Z(i),i=1,ntot)

* create contact map of b-b
       cnatBB=0
       maxdist=0
       do m = 1, nres-4
        do k = m+4, nres
         DIST = ((x(bbat(m))-x(bbat(k)))**2
     $         + (y(bbat(m))-y(bbat(k)))**2
     $         + (z(bbat(m))-z(bbat(k)))**2)**0.5
         if(DIST .LE. NATCUT) then
          cnatBB = cnatBB + 1
          contactmapBB(m,k) = 1
          contactmapBB(k,m) = 1
          distmapBB(m,k)=DIST
          distmapBB(k,m)=DIST
          DIST = sdist*DIST
          if(DIST .GT. maxdist) maxdist=dist
         endif
        enddo
       enddo

* calculate qbb-ns for each structural element
       do i = 1, nstruc 
        cnatsbb(i)=0
        do m = struc(i,1), struc(i,2)
         do k = 1, nres
          if(k > struc(i,2) .or. k < struc(i,1)) then
           if(contactmapBB(m,k) == 1) cnatsbb(i)=cnatsbb(i)+1
          endif
         enddo
        enddo
*        print*,'# of contacts for stuc', i, cnatsbb(i)
       enddo

* calculate number of Type 1 contacts in reference state
       typeonetot=0
       do i = 1, nstruc
        do m = struc(i,1), struc(i,2)
         do k = m+1, struc(i,2)
          if(contactmapBB(m,k) == 1) typeonetot=typeonetot+1
*          print*,'I: ', i, m, k, contactmapBB(m,k), typeonetot
         enddo         
        enddo              
       enddo              

       print*,'# Type I contacts in native state:',typeonetot

* calculate number of Type 2 contacts in reference state
       typetwotot=0                
       do i = 1, nstruc
        do m = struc(i,1), struc(i,2)  
         do f = i+1, nstruc
          do k = struc(f,1), struc(f,2)
           if(contactmapBB(m,k) == 1) typetwotot=typetwotot+1            
*           PRINT*,'II: ',i,m,f,k,contactmapBB(m,k),typetwotot
          enddo 
         enddo
        enddo
       enddo 

       print*,'# Type II contacts in native state:',typetwotot

* calculate number of Type 3 contacts in reference state
       typethreetot=0
       do i = 1, nstruc
        do m = struc(i,1), struc(i,2)
         do k = 1, nres
          use_num = 1
          do f = 1, nstruc
           if (k <= struc(f,2) .and. k >= struc(f,1)) then
            use_num = 0
           endif
          enddo
          if ((contactmapBB(m,k) == 1) .and. use_num == 1) then
           typethreetot=typethreetot+1
          endif
         enddo
        enddo
       enddo

       if(model .eq. 'cacb') then
* create contact map of b-sc
        cnatBS=0
        do m = 1, nres
         do k = 1, nres
          if((scat(k).ne.0).and.((m.lt.k-3).or.(m.gt.k+3))) then
           DIST = ((x(bbat(m))-x(scat(k)))**2
     $         + (y(bbat(m))-y(scat(k)))**2
     $         + (z(bbat(m))-z(scat(k)))**2)**0.5
           if((DIST .LE. NATCUT).and.((m.lt.k-3).or.(m.gt.k+3))) then
            cnatBS = cnatBS + 1
            contactmapBS(m,k) = 1
            distmapBS(m,k)=DIST
            DIST = sdist*DIST
            if(DIST .GT. maxdist) maxdist=dist
           endif
          else
           contactmapBS(m,k) = -9999
          endif
         enddo
        enddo
* create contact map of sc-sc
        cnatSS=0
        do m = 1, nres-3
         do k = m+3, nres
          if((scat(m) .ne. 0) .and. (scat(k) .ne. 0)) then
           DIST = ((x(scat(m))-x(scat(k)))**2
     $         + (y(scat(m))-y(scat(k)))**2
     $         + (z(scat(m))-z(scat(k)))**2)**0.5
           if(DIST .LE. NATCUT) then
            cnatSS = cnatSS + 1
            contactmapSS(m,k) = 1
            distmapSS(m,k)=DIST
            DIST = sdist*DIST
            if(DIST .GT. maxdist) maxdist=dist
           endif
          else
           contactmapSS(m,k) = -9999
          endif
         enddo
        enddo
       endif
 
       cnattot = cnatBB+cnatBS+cnatSS
       PRINT*, cnatBB,cnatBS,cnatSS
      endif

* rg native state
      Rcm(1) = 0
      Rcm(2) = 0
      Rcm(3) = 0
      do m = 1, nat
       Rcm(1) = Rcm(1) + X(m)
       Rcm(2) = Rcm(2) + Y(m)
       Rcm(3) = Rcm(3) + Z(m)
      enddo
      Rcm(1) = Rcm(1)/dble(nat)
      Rcm(2) = Rcm(2)/dble(nat)
      Rcm(3) = Rcm(3)/dble(nat)

      Rg = 0
      do j = 1,nat
       Rg = Rg + (X(J)-Rcm(1))**2
     $ + (Y(J)-Rcm(2))**2 + (Z(J)-Rcm(3))**2
      enddo
      RgNS = SQRT(Rg/dble(nat))

* read in the coordinates
      NFRAMES=ICNTRL8(1)
      PRINT*,'# Frames=',NFRAMES
* do analysis for only the first frame
      NFRAMES = 1  
      DO I=1,NFRAMES
        READ (56) (X(J),J=1,ntot)
        READ (56) (Y(J),J=1,ntot)
        READ (56) (Z(J),J=1,ntot)

* compute number bb native contacts
        cinatBB = 0
        cinnatBB = 0
        do m = 1, nres-4
         do k = m+4, nres
          DIST = ((x(bbat(m))-x(bbat(k)))**2
     $         + (y(bbat(m))-y(bbat(k)))**2
     $         + (z(bbat(m))-z(bbat(k)))**2)**0.5
          if(model .eq. 'cacb') then
           if((DIST .LE. NATCUT).and.(contactmapBB(m,k) /= 1)) then
            cinnatBB = cinnatBB + 1
           endif
           if((DIST.LE.sdist*distmapBB(m,k)).and.
     $        (contactmapBB(m,k) == 1)) then
            cinatBB = cinatBB + 1      
           endif
          else
           if((DIST .LE. NATCUT).and.(contactmapBB(m,k) /= 1)) then
            cinnatBB = cinnatBB + 1
           endif
           if((DIST.LE.sdist*distmapBB(m,k)).and.
     $        (contactmapBB(m,k) == 1)) then
            cinatBB = cinatBB + 1   
           endif  
          endif
         enddo
        enddo
        cinattot=cinatBB
        cinnattot=cinnatBB

        do j = 1, nstruc
         cinatsbb(j)=0
         do m = struc(j,1), struc(j,2)
          do k = 1, nres
           if(k > struc(j,2) .or. k < struc(j,1)) then
            if(contactmapBB(m,k) == 1) then
             DIST = ((x(bbat(m))-x(bbat(k)))**2
     $         + (y(bbat(m))-y(bbat(k)))**2
     $         + (z(bbat(m))-z(bbat(k)))**2)**0.5
             if(DIST.LE.sdist*distmapBB(m,k)) cinatsbb(j)=cinatsbb(j)+1  
            endif
           endif
          enddo
         enddo
        enddo

* calculate number of Type 1 contacts in frame i
        typeonei=0
        do j = 1, nstruc
         do m = struc(j,1), struc(j,2)
          do k = m+1, struc(j,2)
           if(contactmapBB(m,k) == 1) then
            DIST = ((x(bbat(m))-x(bbat(k)))**2
     $       + (y(bbat(m))-y(bbat(k)))**2
     $       + (z(bbat(m))-z(bbat(k)))**2)**0.5
            if (DIST.LE.sdist*distmapBB(m,k)) then
             typeonei=typeonei+1
            endif
           endif
*           print*,'I:',j,m,k,contactmapBB(m,k),typeonei 
          enddo
         enddo
        enddo

* calculate number of Type 2 contacts in frame i
        typetwoi=0
        do j = 1, nstruc
         do m = struc(j,1), struc(j,2)
          do f = j+1, nstruc
           do k = struc(f,1), struc(f,2)
            if(contactmapBB(m,k) == 1) then
             DIST = ((x(bbat(m))-x(bbat(k)))**2 
     $        + (y(bbat(m))-y(bbat(k)))**2  
     $        + (z(bbat(m))-z(bbat(k)))**2)**0.5  
             if (DIST.LE.sdist*distmapBB(m,k)) then
              typetwoi=typetwoi+1
             endif
            endif
*            PRINT*,'II:',j,m,f,k,contactmapBB(m,k),typetwoi
           enddo
          enddo
         enddo
        enddo

* calculate number of Type 3 contacts in frame i
        typethreei=0
        do j = 1, nstruc
         do m = struc(j,1), struc(j,2)
          do k = 1, nres
           use_num = 1
           do f = 1, nstruc
            if (k <= struc(f,2) .and. k >= struc(f,1)) then
             use_num = 0
            endif
*            print*,'Type III', j, m, k, f, use_num
           enddo
           if ((contactmapBB(m,k) == 1) .and. use_num == 1) then
            DIST = ((x(bbat(m))-x(bbat(k)))**2
     $       + (y(bbat(m))-y(bbat(k)))**2
     $       + (z(bbat(m))-z(bbat(k)))**2)**0.5
            if (DIST.LE.sdist*distmapBB(m,k)) then
             typethreei=typethreei+1
            endif
           endif
          enddo
         enddo
        enddo

        if((model .eq. 'cacb').and.(NSREF .ne. 'none')) then              
* compute number bs native contacts
         cinatBS = 0
         cinnatBS = 0
         do m = 1, nres
          do k = 1, nres
           if(contactmapBS(m,k) .NE. -9999) then
            DIST = ((x(bbat(m))-x(scat(k)))**2
     $         + (y(bbat(m))-y(scat(k)))**2
     $         + (z(bbat(m))-z(scat(k)))**2)**0.5
            if((DIST .LE. NATCUT).and.(contactmapBS(m,k) /= 1)) then
             cinnatBS = cinnatBS + 1
            endif
            if((DIST.LE.sdist*distmapBS(m,k)).and.
     $         (contactmapBS(m,k) == 1)) then
              cinatBS = cinatBS + 1
            endif
           endif
          enddo
         enddo
* compute number ss native contacts
         cinatSS = 0
         cinnatSS = 0
         do m = 1, nres-2
          do k = m+2, nres
           if(contactmapSS(m,k) .NE. -9999) then
            DIST = ((x(scat(m))-x(scat(k)))**2
     $         + (y(scat(m))-y(scat(k)))**2
     $         + (z(scat(m))-z(scat(k)))**2)**0.5
            if((DIST .LE. NATCUT).and.(contactmapSS(m,k) /= 1)) then
             cinnatSS = cinnatSS + 1
            endif
            if((DIST.LE.sdist*distmapSS(m,k)).and.
     $         (contactmapSS(m,k) == 1)) then
             cinatSS = cinatSS + 1
            endif
           endif
          enddo
         enddo
         cinattot=cinattot+cinatBS+cinatSS
         cinnattot=cinnattot+cinnatBS+cinnatSS
        endif

* compute rg
        Rcm(1) = 0
        Rcm(2) = 0
        Rcm(3) = 0
        do m = 1, nat
         Rcm(1) = Rcm(1) + X(m)
         Rcm(2) = Rcm(2) + Y(m)
         Rcm(3) = Rcm(3) + Z(m)
        enddo
        Rcm(1) = Rcm(1)/nat
        Rcm(2) = Rcm(2)/nat
        Rcm(3) = Rcm(3)/nat

        Rg = 0
        do j = 1,nat
         Rg = Rg + (X(J)-Rcm(1))**2
     $   + (Y(J)-Rcm(2))**2 + (Z(J)-Rcm(3))**2
        enddo
        Rg = SQRT(Rg/nat)

        e2ed = ((x(bbat(1))-x(bbat(nres)))**2
     $         + (y(bbat(1))-y(bbat(nres)))**2
     $         + (z(bbat(1))-z(bbat(nres)))**2)**0.5

        if((model .eq. 'cacb').and.(NSREF .ne. 'none')) then 

115     FORMAT(I10, 4F10.4, I5, I5,  F8.4, I5, I5, F8.4, I5, I5,  F8.4, 
     $    I5, I5, F8.4, I5, I5)
        WRITE (57,FMT=115) I,Rg/RgNS,Rg,e2ed,
     $   Real(cinatBB)/Real(cnatBB),cinatBB,cinnatBB,
     $   Real(cinatBS)/Real(cnatBS),cinatBS,cinnatBS,
     $   Real(cinatSS)/Real(cnatSS),cinatSS,cinnatSS,
     $   Real(cinattot)/Real(cnattot),cinattot,cinnattot

        elseif((model .eq. 'ca').and.(NSREF .ne. 'none')) then

116     FORMAT(F8.4)
        WRITE (57,FMT=116)
     $   Real(typeonei+typetwoi)/Real(typeonetot+typetwotot)

        elseif((model .eq. 'ca').and.(NSREF .eq. 'none')) then

117     FORMAT(I10, 2F10.4, I5)
        WRITE (57,FMT=117) I,Rg,e2ed,cinnatBB

        endif

* frame read in
      ENDDO
      close(56)
      close(57)

      if((model .eq. 'ca').and.(NSREF .eq. 'none')) then
* write out contact map
* open output data file
      do m = 1, nres
       do k = 1, nres
        if(k<m) then
         pij(m,k)=Real(contactmapBB(k,m))/Real(NFRAMES)
        else
         pij(m,k)=Real(contactmapBB(m,k))/Real(NFRAMES)
        endif
       enddo
      enddo

      OPEN(59,FILE='contact_map.dat', STATUS='UNKNOWN',
     $      FORM='formatted')
       do m = 1, nres
118     FORMAT(200F7.3)
        write(59,118) (pij(m,k),k=1,nres)
       enddo
    
      close(59)
      endif

      END
