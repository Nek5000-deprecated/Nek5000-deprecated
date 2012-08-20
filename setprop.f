      subroutine setprop
C------------------------------------------------------------------------
C
C     Set variable property arrays
C
C------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
C
C     Caution: 2nd and 3rd strainrate invariants residing in scratch
C              common /SCREV/ are used in STNRINV and NEKASGN
C
      common /screv/ sii (lx1,ly1,lz1,lelt),siii(lx1,ly1,lz1,lelt)
      integer offset,ptr

#ifndef NOTIMER
      if (icalld.eq.0) tspro=0.0
      icalld=icalld+1
      nspro=icalld
      etime1=dnekclock()
#endif

      NXYZ1 = NX1*NY1*NZ1
      MFIELD=2
      IF (IFFLOW) MFIELD=1
      nfldt = nfield
      if (ifmhd) nfldt = nfield+1

      ifld = ifield

      DO IFIELD=MFIELD,nfldt
         IF (IFSTRS .AND. IFIELD.EQ.1) CALL STNRINV ! expensive !

         CALL VPROPS

       if (ifadapt) then
           ntot1 =  ntota(IFIELD)
           ptr = adptr(1,ifield)
           if (ifield.ne.1) ptr = ptr+ntota(1)
         avdiff(ifield) = glsc2 (bm1a(ptr),vdiffa(ptr),
     $                           ntot1)/vol
         avtran(ifield) = glsc2 (bm1a(ptr),vtransa(ptr),
     $                           ntot1)/vol
         
        else 
         nel = nelfld(ifield)
         vol = volfld(ifield)
         ntot1 = nxyz1*nel
         avdiff(ifield) = glsc2 (bm1,vdiff (1,1,1,1,ifield),ntot1)/vol
         avtran(ifield) = glsc2 (bm1,vtrans(1,1,1,1,ifield),ntot1)/vol
        endif

      ENDDO

      ifield = ifld

#ifndef NOTIMER
      tspro=tspro+(dnekclock()-etime1)
#endif

C
      RETURN
      END

