      subroutine adinit

      include 'SIZE'
      include 'ADAPT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'GEOM'

      integer nxyz,i,j,lxyz,nel,ptr,ptr1,ptr2

      debug = .false.

      ! Initialize pointers
      nxyz2 = 2*ndim*nx1*ny1 
      nxyz = nx1*ny1*nz1
      lxyz = lx1u*ly1u*lz1u
      do i=1,ldimt1
        nel = nelfld(i)
        ntota(i) = nxyz*nel
        maxord(i) = nx1
        minord(i) = nx1
        do j=1,nel
          adptr(j,i) = 1+nxyz*(j-1)
          bptr(j,i) = 1+nxyz2*(j-1+(i-1)*nel)
        end do
      end do
      call icopy(adptr1(1,1),adptr(1,1),lelt*ldimt1)
      do i=3,ldimt1
        call iadd(adptr(1,i),ntota(1)*(i-2),nelt) 
      end do

!      print *, 'NTOTA'
!      print '(32I5)', ntota
!      print *
!      print *, 'ADPTR'
!      print '(32I5)', adptr
!      print *
!      print *, 'ADPTR1'
!      print '(32I5)', adptr1
!      print *

       nxd = lxdu
       nyd = lydu
       nzd = lzdu

      do j=1,ldimt1
        do i=1,lelt
          nx1var(i,j) = nx1
          ny1var(i,j) = ny1
          nz1var(i,j) = nz1
          nxdvar(i,j) = nxd
          nydvar(i,j) = nyd
          nzdvar(i,j) = nzd
        end do
      end do

      ! Copy initial field values and properties to adaptive solution arrays 
      
      call rzero(vgradt1a,lxyz*lelt*ldimt)
      call rzero(vgradt2a,lxyz*lelt*ldimt)
      ! Mass arrays
       do i=1,ldimt1
          ptr = adptr(1,i)
          nel = nelfld(i)
          if (i.eq.1) then
            call copy(bm1a(ptr),bm1,nxyz*nelv)
            call copy(binvm1a(ptr),binvm1,nxyz*nelv)
            call copy(vxa,vx,nxyz*nelv)
            call copy(vya,vy,nxyz*nelv)
            call copy(vza,vz,nxyz*nelv)
            call copy(vtransa,vtrans,nxyz*nelv)
            call copy(vdiffa,vdiff,nxyz*nelv)
            call copy(v1maska,v1mask,nxyz*nelv)
            call copy(v2maska,v2mask,nxyz*nelv)
            call copy(v3maska,v3mask,nxyz*nelv)
            call copy(vmulta,vmult,nxyz*nelt)
          else
            call copy(Tad(ptr),T(1,1,1,1,i-1),nxyz*nel) 
            call copy(tmaska(ptr),tmask(1,1,1,1,i-1),nxyz*nel)
            call copy(tmulta(ptr),tmult(1,1,1,1,i-1),nxyz*nel)
            do j=1,lorder-1
              ptr2=ivlsum(ntota(2),(i-2))*(lorder-1)
     $             + ntota(i)*(j-1)+1
             call copy(tlaga(ptr2),tlag(1,1,1,1,j,i-1),
     $                 nxyz*nel) 
            end do
            ptr = ptr+ntota(1)
            call copy(bm1a(ptr),bm1,nxyz*nel)
            call copy(binvm1a(ptr),bintm1,nxyz*nel)
            call copy(vtransa(ptr),vtrans(1,1,1,1,i),nxyz*nel)
            call copy(vdiffa(ptr),vdiff(1,1,1,1,i),nxyz*nel)
          end if
          ptr = bptr(1,i)
          call copy(areaa(ptr),area,nxyz2*nel)
!          call copy(xm1a(ptr),xm1,nxyz*nelt)
!          call copy(ym1a(ptr),ym1,nxyz*nelt)
!          call copy(zm1a(ptr),zm1,nxyz*nelt)
       end do
      

      ! Initialize
 
!      call setupds(gsh_fld(ldimt1+1),3,3,3,nelv,
!     $             nelgv,vertex,glo_num)
!
!      call ifill(adords,nx1,27*lelt*ldimt1) 

      return
      end
C----------------------------------------------------------------------
      subroutine padapt

      include 'SIZE'
      include 'ADAPT'
      include 'WZ'
      include 'SOLN'     
      include 'TSTEP'
      include 'PARALLEL' 
      include 'GEOM'
 
      real err(5,2),errtol,wk1(lx1u),wk2(lx1u*ly1u),wk3(lx1u*ly1u*lz1u),
     $     uh(lx1u*ly1u*lz1u),Lj(lx1u*lx1u),Ljt(lx1u*lx1u)
      integer ie,i,ifld,iel,ptr,ptr1,ptr2
      logical if3d
      common /c_is1/ glo_num(1*lx1u*ly1u*lz1u*lelv)
      integer*8 glo_num,ngv
      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex,nxo
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal


      ! Need temporary velocity vector for interpolating and updating
      if3d = .false.
      if (ldim==3) if3d = .true.

      errtol = 1.0D-2
      
      ! Order increment
      adinc = 1 
      
      do ifld=1,ldimt1
        do ie=1,nelt
          adaptflag(ie,ifld) = .false.
        end do
        adaptfld(ifld) = .false.
      end do
      adaptany = .false.
      
      do ifld=2,ldimt1
      nx1 = 0
      do ie=1,nelt
        nxo = nx1
        call getord(ie,ifld)
        ptr = adptr(ie,ifld)
!        if (nx1.ne.lx1u) then
        !  if (nx1.ne.nxo) then
        !    call zwgll(wk1,wk2,nx1)  ! wk1 holds the nodes
        !    call build_legend_transform(Lj,Ljt,wk1,nx1)
        !  end if
          !if (ifld.eq.1) then
          !  call local_err_est(err,vxa(ptr),nx1,Lj,Ljt,uh,wk3,if3d)
          !else 
          !  call local_err_est(err,Tad(ptr),nx1,Lj,Ljt,uh,wk3,if3d)
          !end if

          err(5,1) = 0.0
          if (istep==11.and.(ifld>1)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==41.and.(ifld>1)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==13.and.(ifld>3)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==14.and.(ifld>4)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==15.and.(ifld>5)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==16.and.(ifld>6)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (istep==17.and.(ifld>7)
     $      .and.(nx1+adinc).le.lx1u) err(5,1)=1.0
          if (err(5,1)>errtol) then
             adaptflag(ie,ifld) = .true.           
             adaptfld(ifld) = .true.
             adaptany = .true.
             maxord(ifld) = maxord(ifld)+adinc
             minord(ifld) = minord(ifld)+adinc
          end if 
!        end if
       end do
      end do

      ! Make new pointers
      if (adaptany) then
        newbptr(1,1) = 1
        do ifld=1,ldimt1
          nel = nelfld(ifld)
          newptr1(1,ifld) = 1
          newptr(1,ifld) = 1
  
          if (ifld.gt.2) newptr(1,ifld) = 1+ivlsum(newntota(2),ifld-2)
          do ie=1,nel-1
            call getord(ie,ifld)
            if (adaptflag(ie,ifld)) then
              if (if3d) then
                nxyz = (nx1+adinc)*(ny1+adinc)*(nz1+adinc)
                nxyz2 = (nx1+adinc)*(nz1+adinc)*6
              else
                nxyz = (nx1+adinc)*(ny1+adinc)
                nxyz2 = (nx1+adinc)*4
              end if
            else
              nxyz = nx1*ny1*nz1
              nxyz2 = nx1*nz1*2*ndim
            end if
              newptr(ie+1,ifld) = newptr(ie,ifld)+nxyz
              newptr1(ie+1,ifld) = newptr1(ie,ifld)+nxyz
              newbptr(ie+1,ifld) = newbptr(ie,ifld)+nxyz2
          end do
          if (adaptflag(nel,ifld)) then
            if (if3d) then
              nxyz = (nx1+adinc)*(ny1+adinc)*(nz1+adinc)
              nxyz2 = (nx1+adinc)*(nz1+adinc)*6
            else
              nxyz = (nx1+adinc)*(ny1+adinc)
              nxyz2 = (nx1+adinc)*4
            end if
          else
            nxyz = nx1*ny1*nz1
            nxyz2 = nx1*nz1*2*ndim
          end if
          newntota(ifld) = newptr1(nel,ifld)+nxyz-1
          if (ifld.ne.ldimt1) newbptr(1,ifld+1) =
     $                        newbptr(nel,ifld)+nxyz2
        end do

        !print *, 'PTRS'
        !print '(15I8)', bptr
        !print *
        !print '(15I8)', newbptr
        !print *

        !print *, 'Update Ptrs'
        !print *, 'ptr'
        !print '(6I8)', adptr
        !print '(6I8)', newptr
        !print *, 'ptr1'
        !print '(6I8)', adptr1
        !print '(6I8)', newptr1
        !print *, 'ntot'
        !print '(6I8)', ntota
        !print '(6I8)', newntota
        !print *

        ! Update scalars
        !call UpdateXYZ   ! Do the area update here as well
        call UpdateSol
        !print *, 'AREAA first'
        !print '(16F8.3)', areaa
        !print *
        call UpdateMass
        !print *, 'AREAA second'
        !print '(16F8.3)', areaa
        !print *
   
        ! Multiply vgradt1 and vgradt2 by new mass matrix
        ntot = ivlsum(newntota(2),ldimt)
        ptr = newptr(1,2)
        call col2(vgradt1a,bm1a(ptr+newntota(1)),ntot)
        call col2(vgradt2a,bm1a(ptr+newntota(1)),ntot)
 
        ! Update number of nodes 
        do ifld=1,ldimt1
          nel = nelfld(ifld) 
          do ie=1,lelt
            if (adaptflag(ie,ifld)) then
              nx1var(ie,ifld) = nx1var(ie,ifld)+adinc
              ny1var(ie,ifld) = ny1var(ie,ifld)+adinc
              !nxdvar(ie,ifld) = 3*nx1var(ie,ifld)/2
              !nydvar(ie,ifld) = 3*ny1var(ie,ifld)/2
              if (if3d) then
                nz1var(ie,ifld) = nz1var(ie,ifld)+adinc
              !  nzdvar(ie,ifld) = 3*nz1var(ie,ifld)/2
              end if
              !call ifill(adords(1,1,1,ie,ifld),nx1,27)
            end if
          end do
        end do
  
        !print *, 'ADPTR NEWPTR', nid
        !print *, adptr
        !print *
        !print *, newptr
        !print *
       
        !print *, 'ADPTR1 NEWPTR1', nid
        !print *, adptr1
        !print *
        !print *, newptr1
        !print *
       
        !print *, 'NTOTA', nid
        !print *, ntota
        !print *
        !print *, newntota
        !print * 
 
        call icopy(adptr,newptr,lelt*ldimt1)
        call icopy(adptr1,newptr1,lelt*ldimt1)
        call icopy(ntota,newntota,ldimt1)
        call icopy(bptr,newbptr,lelt*ldimt1)

        ! Update DS
        do ifld=1,ldimt1
          if (adaptfld(ifld)) then  ! Need to do an all to all so all procs have the correct flag. 
            nel = nelfld(ifld)
            call izero(glo_num,lx1u*ly1u*lz1u*nel)
            call getord(1,ifld)
       call setupds(gsh_fld(ifld),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
!            call sv3d_modal(glo_num,ngv,nx1var,ny1var,nz1var,nel,
!     $                      vertex,.false.)
!
!            call gs_setup(gsh_fld(ifld),glo_num,ntota(ifld),nekcomm,mp)
            !print *, 'GLO_NUM'
            !if (nx1==4) then
            !  print '(16I6)', glo_num
            !elseif (nx1==5) then 
            !  print '(25I6)', glo_num
            !end if
            !print *
            !print *, 'GLO_NUM', nx1
            !if (nx1.eq.3) then
            !  print '(9I6)', glo_num
            !elseif (nx1.eq.4) then
            !  print '(16I6)', glo_num
            !else
            !  print '(25I6)', glo_num
            !end if
          end if
        end do

        ! Update Mult
         ! print *, 'TMULT', nx1
         ! if (nx1.eq.3) then
         !   print '(25F8.3)', tmulta 
         ! elseif (nx1.eq.4) then
         !   print '(25F8.3)', tmulta
         ! else
         !   print '(25F8.3)', tmulta
         ! end if
         ! print *
        call UpdateMult
        ! Update Mask
         ! if (nx1.eq.3) then
         !   print '(25F8.3)', tmulta 
         ! elseif (nx1.eq.4) then
         !   print '(25F8.3)', tmulta
         ! else
         !   print '(25F8.3)', tmulta
         ! end if
         ! print *
         ! print *, 'TMASK'
         ! print '(49F8.3)', tmaska
         ! print *
        call UpdateBCs
         ! print *, 'TMASK 2'
         ! print '(64F8.3)', tmaska
         ! print *
        !print *
        !print *, 'PROP 1'
        !print '(16F8.3)', vdiffa
        !print *
        !print '(16F8.3)', vtransa
        !print *
        call setprop
        !print *
        !print *, 'PROP 2'
        !print '(25F8.3)', vdiffa
        !print *
        !print '(25F8.3)', vtransa
        !print *
      end if

      return
      end
c-----------------------------------------------------------------------
      subroutine getord(ie,ifld)
  
      include 'SIZE' 
      include 'ADAPT'

      integer ie,ifld

      nx1 = nx1var(ie,ifld)
      ny1 = ny1var(ie,ifld)
      nz1 = nz1var(ie,ifld)
      nxd = (2*nx1+lx1)/2 
      nyd = (2*ny1+ly1)/2 
      nzd = (2*nz1+lz1)/2 
!      nxd = nxdvar(ie,ifld)
!      nyd = nydvar(ie,ifld)
!      nzd = nzdvar(ie,ifld)
   
      end
c----------------------------------------------------------------------
      subroutine setmapa(nxold,nxnew)

      include 'SIZE'
      include 'ADAPT'
      include 'WZ'

      integer nxnew, nxold

      ! Get new nodes, derivative matrices, and interpolation matrices
      real zgm1o(lx1u*3),zgm1a(lx1u*3)
 
!      if (ndim.eq.2) then
        ! Old nodes
        CALL ZWGLL (ZGM1O(1),WXM1E,NXOLD)
        !CALL ZWGLL (ZGM1O(nx1+1),WYM1E,NY1)
        ! New nodes
        CALL ZWGLL (ZGM1A(1),WXM1E,NXNEW)
        !CALL ZWGLL (ZGM1A(nx1+1),WYM1E,NYNEW)
        CALL IGLLM(IXADPT,IXTADPT,zgm1o(1),zgm1a(1),NXOLD,NXNEW,
     $             NXOLD,NXNEW)
        !CALL IGLLM(IYADPT,IYTADPT,zgm1o(nx1+1),zgm1a(nx1+2),ny1,ny1+1,
     $  !           ny1,ny1+1)
        
!      else
       ! CALL ZWGLL (ZGM1A(1,1),WXM1E,NX1)
       ! CALL ZWGLL (ZGM1A(1,2),WYM1E,NY1)
       ! CALL ZWGLL (ZGM1A(1,3),WZM1E,NZ1)
       ! CALL IGLLM(IXADPT,IXTADPT,zgm1o(1,1),zgm1a(1,1),nx1-1,nx1,
     $ !            nx1-1,nx1)
       ! CALL IGLLM(IYADPT,IYTADPT,zgm1o(1,2),zgm1a(1,2),ny1-1,ny1,
     $ !            ny1-1,ny1)
       ! CALL IGLLM(IZADPT,IZTADPT,zgm1o(1,3),zgm1a(1,3),nz1-1,nz1,
     $ !            nz1-1,nz1)
!      end if

      end 
c-----------------------------------------------------------------------
      subroutine interp(unew,uold,ifld)

      include "SIZE"
      include "ADAPT"
      include 'INPUT'
      include 'TSTEP'
    
      real uold(1),wk1(lx1u*ly1u*lz1u),
     $     wk2(lx1u*ly1u*lz1u),unew(1)
      integer ie,iord,ntot,ifld,nzold,nxyz,ptr,ptr1,ptr2

      iord = 0 
      nel = nelfld(ifld)
      do ie=1,nel
        ptr1 = adptr1(ie,ifld)
        ptr2 = newptr1(ie,ifld)
        call getord(ie,ifld) 
        if (adaptflag(ie,ifld)) then
          if (nx1.ne.iord) then
            call setmapa(nx1,nx1+adinc) 
            iord = nx1
          end if
          ! Assumes nx1 = ny1 = nz1
          call specmp(unew(ptr2),nx1+adinc,uold(ptr1),
     $               nx1,ixadpt,ixtadpt,wk1)
!          if (if3d) then
!            CALL MXM (ixadpt,NX1+1,uold(ptr1),NX1,wk1,ny1*nz1)
!            DO 5 IZ=1,nz
!            ptr3 = (IZ-1)*(NX1+1)*NY1+1
!            CALL MXM (wk1(ptr3),NX1+1,iytadpt,NY1,wk2,NY1+1)
! 5          CONTINUE
!            CALL MXM (wk2,(NX1-1)*(NY1-1),iztadpt,nz1,unew(ptr2),NZ1+1)
!          else
!            CALL MXM (ixadpt,NX1+1,uold(ptr1),NX1,wk1,ny1)
!            CALL MXM (wk1,NX1+1,iytadpt,ny1,unew(ptr2),ny1+1)
!          end if
        else
          nxyz = nx1*ny1*nz1
          call copy(unew(ptr2),uold(ptr1),nxyz)
        end if
      end do
   
      return 
      end 
c----------------------------------------------------------------------
      subroutine ComputeG(ie,Gvals)
      
      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'TOPOL'
      include 'INPUT'
      include 'PARALLEL'
      include 'ADAPT'    
  
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      real Gvals(nx1*ny1*nz1*6)
      integer ie,nyz1,ptrvec(6),ctr,iel,ptr,ptr1,ptr2
      real  XRM1E(NX1*NY1*NZ1),YRM1E(NX1*NY1*NZ1),XSM1E(NX1*NY1*NZ1),
     $      YSM1E(NX1*NY1*NZ1),XTM1E(NX1*NY1*NZ1),YTM1E(NX1*NY1*NZ1),
     $      ZRM1E(NX1*NY1*NZ1),ZSM1E(NX1*NY1*NZ1),ZTM1E(NX1*NY1*NZ1)
      ! Get x, y, z for this element
      call genxyz1(xm1e,ym1e,zm1e,nx1,ny1,nz1,ie)

      NXY1  = NX1*NY1
      NYZ1  = NY1*NZ1
      NXYZ1 = NX1*NY1*NZ1
      ntot1 = nxyz1
      ! Get new derivative matrices
      call ComputeDxyz 
   
      ! Get partial derivatives for this element
      CALL XYZRST1(XRM1E,YRM1E,ZRM1E,XSM1E,YSM1E,ZSM1E,XTM1E,YTM1E,
     $             ZTM1E,IFAXIS,IE)
    
      do i=1,6
        ptrvec(i) = 1+(i-1)*nxyz1 
      end do

      IF (NDIM.EQ.2) THEN
         CALL RZERO   (JACM1E,NTOT1)
         CALL ADDCOL3 (JACM1E,XRM1E,YSM1E,NTOT1)
         CALL SUBCOL3 (JACM1E,XSM1E,YRM1E,NTOT1)
         CALL COPY    (RXM1E,YSM1E,NTOT1)
         CALL COPY    (RYM1E,XSM1E,NTOT1)
         CALL CHSIGN  (RYM1E,NTOT1)
         CALL COPY    (SXM1E,YRM1E,NTOT1)
         CALL CHSIGN  (SXM1E,NTOT1)
         CALL COPY    (SYM1E,XRM1E,NTOT1)
         CALL RZERO   (RZM1E,NTOT1)
         CALL RZERO   (SZM1E,NTOT1)
         CALL RONE    (TZM1E,NTOT1)
         call invers2(jacmie,jacm1e,ntot1)
c        Compute geometric factors for integrated del-squared operator
         CALL VDOT2 (Gvals(ptrvec(1)),RXM1E,RYM1E,RXM1E,RYM1E,NTOT1)
         CALL VDOT2 (Gvals(ptrvec(2)),SXM1E,SYM1E,SXM1E,SYM1E,NTOT1)
         CALL VDOT2 (Gvals(ptrvec(4)),RXM1E,RYM1E,SXM1E,SYM1E,NTOT1)
         CALL COL2  (Gvals(ptrvec(1)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(2)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(4)),JACMIE,NTOT1)
         CALL RZERO (Gvals(ptrvec(3)),NTOT1)
         CALL RZERO (Gvals(ptrvec(5)),NTOT1)
         CALL RZERO (Gvals(ptrvec(6)),NTOT1)
      ELSE
         CALL RZERO   (JACM1E,NTOT1)
         CALL ADDCOL4 (JACM1E,XRM1E,YSM1E,ZTM1E,NTOT1)
         CALL ADDCOL4 (JACM1E,XTM1E,YRM1E,ZSM1E,NTOT1)
         CALL ADDCOL4 (JACM1E,XSM1E,YTM1E,ZRM1E,NTOT1)
         CALL SUBCOL4 (JACM1E,XRM1E,YTM1E,ZSM1E,NTOT1)
         CALL SUBCOL4 (JACM1E,XSM1E,YRM1E,ZTM1E,NTOT1)
         CALL SUBCOL4 (JACM1E,XTM1E,YSM1E,ZRM1E,NTOT1)
         CALL ASCOL5  (RXM1E,YSM1E,ZTM1E,YTM1E,ZSM1E,NTOT1)
         CALL ASCOL5  (RYM1E,XTM1E,ZSM1E,XSM1E,ZTM1E,NTOT1)
         CALL ASCOL5  (RZM1E,XSM1E,YTM1E,XTM1E,YSM1E,NTOT1)
         CALL ASCOL5  (SXM1E,YTM1E,ZRM1E,YRM1E,ZTM1E,NTOT1)
         CALL ASCOL5  (SYM1E,XRM1E,ZTM1E,XTM1E,ZRM1E,NTOT1)
         CALL ASCOL5  (SZM1E,XTM1E,YRM1E,XRM1E,YTM1E,NTOT1)
         CALL ASCOL5  (TXM1E,YRM1E,ZSM1E,YSM1E,ZRM1E,NTOT1)
         CALL ASCOL5  (TYM1E,XSM1E,ZRM1E,XRM1E,ZSM1E,NTOT1)
         CALL ASCOL5  (TZM1E,XRM1E,YSM1E,XSM1E,YRM1E,NTOT1)
         call invers2 (jacmie,jacm1e,ntot1)
c        Compute geometric factors for integrated del-squared operator
         CALL VDOT3 (Gvals(ptrvec(1)),RXM1E,RYM1E,RZM1E,RXM1E,RYM1E,
     $               RZM1E,NTOT1)
         CALL VDOT3 (Gvals(ptrvec(2)),SXM1E,SYM1E,SZM1E,SXM1E,SYM1E,
     $               SZM1E,NTOT1)
         CALL VDOT3 (Gvals(ptrvec(3)),TXM1E,TYM1E,TZM1E,TXM1E,TYM1E,
     $               TZM1E,NTOT1)
         CALL VDOT3 (Gvals(ptrvec(4)),RXM1E,RYM1E,RZM1E,SXM1E,SYM1E,
     $               SZM1E,NTOT1)
         CALL VDOT3 (Gvals(ptrvec(5)),RXM1E,RYM1E,RZM1E,TXM1E,TYM1E,
     $               TZM1E,NTOT1)
         CALL VDOT3 (Gvals(ptrvec(6)),SXM1E,SYM1E,SZM1E,TXM1E,TYM1E,
     $               TZM1E,NTOT1)
         CALL COL2  (Gvals(ptrvec(1)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(2)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(3)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(4)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(5)),JACMIE,NTOT1)
         CALL COL2  (Gvals(ptrvec(6)),JACMIE,NTOT1)
      ENDIF

C     Multiply the geometric factors GiM1,i=1,5 with the
C     weights on mesh M1.
      IF (IFAXIS) CALL SETAXW1 ( IFRZER(IE) )
      CALL COL2 (Gvals(ptrvec(1)),W3M1E,NXYZ1)
      CALL COL2 (Gvals(ptrvec(2)),W3M1E,NXYZ1)
      CALL COL2 (Gvals(ptrvec(4)),W3M1E,NXYZ1)
      IF (NDIM.EQ.3) THEN
         CALL COL2 (Gvals(ptrvec(3)),W3M1E,NXYZ1)
         CALL COL2 (Gvals(ptrvec(5)),W3M1E,NXYZ1)
         CALL COL2 (Gvals(ptrvec(6)),W3M1E,NXYZ1)
      ENDIF

      IF (NDIM.EQ.3) THEN
            IF (.NOT.IFDFRM(IE)) THEN
               ctr = 0
               DO 1000 IZ=1,NZ1
               DO 1000 IY=1,NY1
               DO 1000 IX=1,NX1
                  Gvals(ptrvec(4)+ctr)=Gvals(ptrvec(1)+ctr)/WXM1E(IX)
                  Gvals(ptrvec(5)+ctr)=Gvals(ptrvec(2)+ctr)/WYM1E(IY)
                  Gvals(ptrvec(6)+ctr)=Gvals(ptrvec(3)+ctr)/WZM1E(IZ)
                  ctr = ctr + 1
 1000          CONTINUE
            ENDIF
      ELSE
            IF (.NOT.IFDFRM(IE)) THEN
               ctr = 0
               DO 2000 IY=1,NY1
               DO 2000 IX=1,NX1                
                  Gvals(ptrvec(4)+ctr)=Gvals(ptrvec(1)+ctr)/WXM1E(IX)
                  Gvals(ptrvec(5)+ctr)=Gvals(ptrvec(2)+ctr)/WYM1E(IY)
                  ctr = ctr+1
 2000          CONTINUE
            ENDIF
      ENDIF

      RETURN
      end
c-----------------------------------------------------------------------
      subroutine genxyz1 (xml,yml,zml,nxl,nyl,nzl,e)
C
      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'TOPOL'
      include 'INPUT'
      include 'PARALLEL'
      include 'ADAPT'

      real xml(nxl*nyl*nzl),yml(nxl*nyl*nzl),zml(nxl*nyl*nzl)

C     Note : CTMP1 is used in this format in several subsequent routines
      common /ctmp1/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1),
     $               zgml(lx1,3),work(3,lx1,lz1)

      parameter (ldw=2*lx1*ly1*lz1)
      common /ctmp0/ w(ldw)

      character*1 ccv
      integer e, nxl, nyl, nzl

#ifdef MOAB
c already read/initialized vertex positions
      if (ifmoab) return
#endif
        
c     Initialize geometry arrays with bi- triquadratic deformations
      call linquad1(xml,yml,zml,nxl,nyl,nzl,e)

      call setzgml (zgml,e,nxl,nyl,nzl,ifaxis)
      call sethmat (h,zgml,nxl,nyl,nzl)

c     Deform surfaces - general 3D deformations
c                    - extruded geometry deformations
      nfaces = 2*ndim
      do iface=1,nfaces
        ccv = ccurve(iface,e)
        if (ccv.eq.'s') 
     $     call sphsrf(xml,yml,zml,iface,e,nxl,nyl,nzl,work) 
        if (ccv.eq.'e') 
     $     call gensrf(xml,yml,zml,iface,e,nxl,nyl,nzl,zgml) 
      enddo

      do isid=1,8
        ccv = ccurve(isid,e)
        if (ccv.eq.'C') then 
        call arcsrf1(xml,yml,zml,nxl,nyl,nzl,e,isid)
        end if
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine linquad1(xl,yl,zl,nxl,nyl,nzl,e)

      ! Same as linquad but only for a single element

      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'TOPOL'
      include 'INPUT'
      include 'PARALLEL'
      include 'ADAPT'

      real xl(nxl*nyl*nzl),yl(nxl*nyl*nzl),zl(nxl*nyl*nzl)

      integer e,nxl,nyl,nzl
      logical ifmid

      nedge = 4 + 8*(ndim-2)

      ifmid = .false.
      do k=1,nedge
         if (ccurve(k,e).eq.'m') ifmid = .true.
      enddo

      if (lx1.eq.2) ifmid = .false.
      if (ifmid) then
         call xyzquad(xl,yl,zl,nxl,nyl,nzl,e)
      else
         call xyzlin1(xl,yl,zl,nxl,nyl,nzl,e,ifaxis)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine xyzrst1(xrm1e,yrm1e,zrm1e,xsm1e,ysm1e,zsm1e,
     $                   XTM1e,YTM1e,ZTM1e,IFAXIS,ie)
C-----------------------------------------------------------------------
C
C     Compute global-to-local derivatives on mesh 1.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ADAPT'
      INCLUDE 'GEOM'
C
      REAL     XRM1E(LX1U*LY1U*LZ1U),YRM1E(LX1U*LY1U*LZ1U),
     $         ZRM1E(LX1U*LY1U*LZ1U),XSM1E(LX1U*LY1U*LZ1U),
     $         YSM1E(LX1U*LY1U*LZ1U),ZSM1E(LX1U*LY1U*LZ1U),
     $         XTM1E(LX1U*LY1U*LZ1U),YTM1E(LX1U*LY1U*LZ1U),
     $         ZTM1E(LX1U*LY1U*LZ1U)
      LOGICAL IFAXIS
      INTEGER IE,ptr,ptr1,ptr2
C
      NXY1=NX1*NY1
      NYZ1=NY1*NZ1

C
      IF (IFAXIS) CALL SETAXDY ( IFRZER(IE) )
C
      CALL MXM (DXM1E,NX1,XM1E,NX1,XRM1E,NYZ1)
      CALL MXM (DXM1E,NX1,YM1E,NX1,YRM1E,NYZ1)
      CALL MXM (DXM1E,NX1,ZM1E,NX1,ZRM1E,NYZ1)
C
      DO 10 IZ=1,NZ1
      ptr = (IZ-1)*NXY1+1
      CALL MXM (XM1E(ptr),NX1,DYTM1E,NY1,XSM1E(ptr),NY1)
      CALL MXM (YM1E(ptr),NX1,DYTM1E,NY1,YSM1E(ptr),NY1)
      CALL MXM (ZM1E(ptr),NX1,DYTM1E,NY1,ZSM1E(ptr),NY1)
   10 CONTINUE
C
      IF (NDIM.EQ.3) THEN
         CALL MXM (XM1E,NXY1,DZTM1E,NZ1,XTM1E,NZ1)
         CALL MXM (YM1E,NXY1,DZTM1E,NZ1,YTM1E,NZ1)
         CALL MXM (ZM1E,NXY1,DZTM1E,NZ1,ZTM1E,NZ1)
      ELSE
         CALL RZERO (XTM1E,NXY1)
         CALL RZERO (YTM1E,NXY1)
         CALL RONE  (ZTM1E,NXY1)
      ENDIF
C
      RETURN
      END
C------------------------------------------------------------------
      subroutine xyzlin1(xl,yl,zl,nxl,nyl,nzl,e,ifaxl)
c     Generate bi- or trilinear mesh

      include 'SIZE'
      include 'INPUT'
      include 'ADAPT'

      real xl(nxl*nyl*nzl),yl(nxl*nyl*nzl),zl(nxl*nyl*nzl)
      integer e
      logical ifaxl ! local ifaxis specification

c   Preprocessor Corner notation:      Symmetric Corner notation:
c
c           4+-----+3    ^ s                    3+-----+4    ^ s
c           /     /|     |                      /     /|     |
c          /     / |     |                     /     / |     |
c        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
c         |     | /     /                     |     | /     /
c         |     |/     /                      |     |/     /
c        5+-----+6    t                      5+-----+6    t

      integer indx(8)
      save    indx
      data    indx / 1,2,4,3,5,6,8,7 /

      parameter (ldw=4*lx1u*ly1u*lz1u)
      common /ctmp0/ xcb(2,2,2),ycb(2,2,2),zcb(2,2,2),w(ldw)

      common /cxyzla/ zgmla(lx1u*3),jxa(lx1u*2),jya(lx1u*2),jza(lx1u*2)
     $                          ,jxta(lx1u*2),jyta(lx1u*2),jzta(lx1u*2)
     $                          ,zlina(2)
      real jxa,jya,jza,jxta,jyta,jzta

      !call setzgml (zgmla,e,nxl,nyl,nzl,ifaxl)
      call rzero(zgmla,3*lx1u)
      !call zwgll(zgmla(1),nxl,w)
      !call zwgll(zgmla(nxl+1),nyl,w)
      !if (nzl.gt.1) then
      !  call zwgll(zgmla(nxl+nyl+1),nzl,w)
      !end if

      call zwgljd(zgmla(1),w,nxl,0.0,0.0)
      call zwgljd(zgmla(nxl+1),w,nyl,0.0,0.0)
      if (nzl.gt.1) then
        call zwgljd(zgmla(nxl+nyl+1),w,nzl,0.0,0.0)
      end if
     
      zlina(1) = -1
      zlina(2) =  1

      k = 1
      do i=1,nxl
         call fd_weights_full(zgmla(i),zlina,1,0,jxta(k))
         call fd_weights_full(zgmla(nxl+i),zlina,1,0,jyta(k))
         call fd_weights_full(zgmla(2*nxl+i),zlina,1,0,jzta(k))
         k=k+2
      enddo
      call transpose(jxa,nxl,jxta,2)

      ndim2 = 2**ndim
      do ix=1,ndim2          ! Convert prex notation to lexicographical
         i=indx(ix)
         xcb(ix,1,1)=xc(i,e)
         ycb(ix,1,1)=yc(i,e)
         zcb(ix,1,1)=zc(i,e)
      enddo
      
c     Map R-S-T space into physical X-Y-Z space.

      ! NOTE:  Assumes nxl=nyl=nzl !
      call tensr3(xl,nxl,xcb,2,jxa,jyta,jzta,w)
      call tensr3(yl,nxl,ycb,2,jxa,jyta,jzta,w)
      call tensr3(zl,nxl,zcb,2,jxa,jyta,jzta,w)
     
      return
      end
c--------------------------------------------------------------------------------
      subroutine arcsrf1(xml,yml,zml,nxl,nyl,nzl,ie,isid)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'WZ'
C
C     ....note..... CTMP1 is used in this format in several subsequent routines
C
      COMMON /CTMP1/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      DIMENSION XML(NXL,NYL,NZL,1),YML(NXL,NYL,NZL,1),ZML(NXL,NYL,NZL,1)
      LOGICAL IFGLJ
C
      IFGLJ = .FALSE.
      IF (IFAXIS .AND. IFRZER(IE) .AND. (ISID.EQ.2 .OR. ISID.EQ.4)) 
     $IFGLJ = .TRUE.
C
      PT1X  = XC(ISID,IE)
      PT1Y  = YC(ISID,IE)
      IF(ISID.EQ.4) THEN
         PT2X = XC(1,IE)
         PT2Y = YC(1,IE)
      ELSE IF(ISID.EQ.8) THEN
         PT2X = XC(5,IE)
         PT2Y = YC(5,IE)
      ELSE
         PT2X = XC(ISID+1,IE)
         PT2Y = YC(ISID+1,IE)
      ENDIF
C
C     Find slope of perpendicular
      RADIUS=CURVE(1,ISID,IE)
      GAP=SQRT( (PT1X-PT2X)**2 + (PT1Y-PT2Y)**2 )
      IF (ABS(2.0*RADIUS).LE.GAP*1.00001) THEN
         write(6,10) RADIUS,ISID,IE,GAP
   10    FORMAT(//,2X,'ERROR: Too small a radius (',G11.3
     $  ,') specified for side',I2,' of element',I4,':  '
     $  ,G11.3,/,2X,'ABORTING during mesh generation.')
         call exitt
      ENDIF
      XS = PT2Y-PT1Y
      YS = PT1X-PT2X
C     Make length Radius
      XYS=SQRT(XS**2+YS**2)
C     Find Center
      DTHETA = ABS(ASIN(0.5*GAP/RADIUS))
      PT12X  = (PT1X + PT2X)/2.0
      PT12Y  = (PT1Y + PT2Y)/2.0
      XCENN  = PT12X - XS/XYS * RADIUS*COS(DTHETA)
      YCENN  = PT12Y - YS/XYS * RADIUS*COS(DTHETA)
      THETA0 = ATAN2((PT12Y-YCENN),(PT12X-XCENN))
      IF (IFGLJ) THEN
         FAC    = SIGN(1.0,RADIUS)
         THETA1 = THETA0 - FAC*DTHETA
         THETA2 = THETA0 + FAC*DTHETA
      ENDIF
C     Compute perturbation of geometry
      ISID1 = MOD1(ISID,4)
      IF (IFGLJ) THEN
         I1 = ISID/2
         I2 = 2 - ISID/4
         DO 15 IY=1,NYL
           ANG  = H(IY,2,I1)*THETA1 + H(IY,2,I2)*THETA2
           XCRVED(IY)=XCENN + ABS(RADIUS)*COS(ANG)
     $                      - (H(IY,2,I1)*PT1X + H(IY,2,I2)*PT2X)
           YCRVED(IY)=YCENN + ABS(RADIUS) * SIN(ANG)
     $                      - (H(IY,2,I1)*PT1Y + H(IY,2,I2)*PT2Y)
   15    CONTINUE
      ELSE
         DO 20 IX=1,NXL
            IXT=IX
            IF (ISID1.GT.2) IXT=NXL+1-IX
            R=ZGML(IX,1)
            IF (RADIUS.LT.0.0) R=-R
            XCRVED(IXT) = XCENN + ABS(RADIUS) * COS(THETA0 + R*DTHETA)
     $                          - ( H(IX,1,1)*PT1X + H(IX,1,2)*PT2X )
            YCRVED(IXT) = YCENN + ABS(RADIUS) * SIN(THETA0 + R*DTHETA)
     $                          - ( H(IX,1,1)*PT1Y + H(IX,1,2)*PT2Y )
   20    CONTINUE
      ENDIF
C     Points all set, add perturbation to current mesh.
      ISID1 = MOD1(ISID,4)
      ISID1 = EFACE1(ISID1)
      IZT = (ISID-1)/4+1
      IYT = ISID1-2
      IXT = ISID1
      IF (ISID1.LE.2) THEN
         CALL ADDTNSR(XML(1,1,1,1),H(1,1,IXT),XCRVED,H(1,3,IZT)
     $               ,NXL,NYL,NZL)
         CALL ADDTNSR(YML(1,1,1,1),H(1,1,IXT),YCRVED,H(1,3,IZT)
     $               ,NXL,NYL,NZL)
      ELSE
         CALL ADDTNSR(XML(1,1,1,1),XCRVED,H(1,2,IYT),H(1,3,IZT)
     $               ,NXL,NYL,NZL)
         CALL ADDTNSR(YML(1,1,1,1),YCRVED,H(1,2,IYT),H(1,3,IZT)
     $               ,NXL,NYL,NZL)
      ENDIF
      return
      end
c------------------------------------------------------------------
!      SUBROUTINE CUMAX (V1,V2,V3,UMAX)
!C
!      INCLUDE 'SIZE'
!      INCLUDE 'WZ'
!      INCLUDE 'GEOM'
!      INCLUDE 'INPUT'
!C
!      COMMON /SCRNS/ XRM1 (LX1,LY1,LZ1,LELV)
!     $ ,             XSM1 (LX1,LY1,LZ1,LELV)
!     $ ,             XTM1 (LX1,LY1,LZ1,LELV)
!     $ ,             YRM1 (LX1,LY1,LZ1,LELV)
!     $ ,             YSM1 (LX1,LY1,LZ1,LELV)
!     $ ,             YTM1 (LX1,LY1,LZ1,LELV)
!      COMMON /SCRMG/ ZRM1 (LX1,LY1,LZ1,LELV)
!     $ ,             ZSM1 (LX1,LY1,LZ1,LELV)
!     $ ,             ZTM1 (LX1,LY1,LZ1,LELV)
!      COMMON /CTMP1/ U    (LX1,LY1,LZ1,LELV)
!     $ ,             V    (LX1,LY1,LZ1,LELV)
!     $ ,             W    (LX1,LY1,LZ1,LELV)
!      COMMON /CTMP0/ X    (LX1,LY1,LZ1,LELV)
!     $ ,             R    (LX1,LY1,LZ1,LELV)
!      COMMON /DELRST/ DRST(LX1),DRSTI(LX1)
!C
!      DIMENSION V1(LX1,LY1,LZ1,1)
!     $        , V2(LX1,LY1,LZ1,1)
!     $        , V3(LX1,LY1,LZ1,1)
!      DIMENSION U3(3)
!      INTEGER ICALLD
!      SAVE    ICALLD
!      DATA    ICALLD /0/
!C
!      NTOT  = NX1*NY1*NZ1*NELV
!      NTOTL = LX1*LY1*LZ1*LELV
!      NTOTD = NTOTL*NDIM
!C
!C     Compute isoparametric partials.
!C
!      CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,
!     $             IFAXIS)
!C
!C     Compute maximum U/DX
!C
!      IF (ICALLD.EQ.0) THEN
!         ICALLD=1
!         DRST (1)=ABS(ZGM1(2,1)-ZGM1(1,1))
!         DRSTI(1)=1.0/DRST(1)
!         DO 400 I=2,NX1-1
!            DRST (I)=ABS(ZGM1(I+1,1)-ZGM1(I-1,1))/2.0
!            DRSTI(I)=1.0/DRST(I)
! 400     CONTINUE
!         DRST (NX1)=DRST(1)
!         DRSTI(NX1)=1.0/DRST(NX1)
!      ENDIF
!C
!C     Zero out scratch arrays U,V,W for ALL declared elements...
!C
!      CALL RZERO3 (U,V,W,NTOTL)
!C
!      IF (NDIM.EQ.2) THEN
!
!      CALL VDOT2  (U,V1  ,V2  ,RXM1,RYM1,NTOT)
!      CALL VDOT2  (R,RXM1,RYM1,RXM1,RYM1,NTOT)
!      CALL VDOT2  (X,XRM1,YRM1,XRM1,YRM1,NTOT)
!      CALL COL2   (R,X,NTOT)
!      CALL VSQRT  (R,NTOT)
!      CALL INVCOL2(U,R,NTOT)
!C
!      CALL VDOT2  (V,V1  ,V2  ,SXM1,SYM1,NTOT)
!      CALL VDOT2  (R,SXM1,SYM1,SXM1,SYM1,NTOT)
!      CALL VDOT2  (X,XSM1,YSM1,XSM1,YSM1,NTOT)
!      CALL COL2   (R,X,NTOT)
!      CALL VSQRT  (R,NTOT)
!      CALL INVCOL2(V,R,NTOT)
!C
!      ELSE
!C
!      CALL VDOT3  (U,V1  ,V2  ,V3  ,RXM1,RYM1,RZM1,NTOT)
!      CALL VDOT3  (R,RXM1,RYM1,RZM1,RXM1,RYM1,RZM1,NTOT)
!      CALL VDOT3  (X,XRM1,YRM1,ZRM1,XRM1,YRM1,ZRM1,NTOT)
!      CALL COL2   (R,X,NTOT)
!      CALL VSQRT  (R,NTOT)
!      CALL INVCOL2(U,R,NTOT)
!C
!      CALL VDOT3  (V,V1  ,V2  ,V3  ,SXM1,SYM1,SZM1,NTOT)
!      CALL VDOT3  (R,SXM1,SYM1,SZM1,SXM1,SYM1,SZM1,NTOT)
!      CALL VDOT3  (X,XSM1,YSM1,ZSM1,XSM1,YSM1,ZSM1,NTOT)
!      CALL COL2   (R,X,NTOT)
!      CALL VSQRT  (R,NTOT)
!      CALL INVCOL2(V,R,NTOT)
!C
!      CALL VDOT3  (W,V1  ,V2  ,V3  ,TXM1,TYM1,TZM1,NTOT)
!      CALL VDOT3  (R,TXM1,TYM1,TZM1,TXM1,TYM1,TZM1,NTOT)
!      CALL VDOT3  (X,XTM1,YTM1,ZTM1,XTM1,YTM1,ZTM1,NTOT)
!      CALL COL2   (R,X,NTOT)
!      CALL VSQRT  (R,NTOT)
!      CALL INVCOL2(W,R,NTOT)
!C
!      ENDIF
!C
!      DO 500 IE=1,NELV
!      DO 500 IX=1,NX1
!      DO 500 IY=1,NY1
!      DO 500 IZ=1,NZ1
!            U(IX,IY,IZ,IE)=ABS( U(IX,IY,IZ,IE)*DRSTI(IX) )
!            V(IX,IY,IZ,IE)=ABS( V(IX,IY,IZ,IE)*DRSTI(IY) )
!            W(IX,IY,IZ,IE)=ABS( W(IX,IY,IZ,IE)*DRSTI(IZ) )
!  500    CONTINUE
!C
!      U3(1)   = VLMAX(U,NTOT)
!      U3(2)   = VLMAX(V,NTOT)
!      U3(3)   = VLMAX(W,NTOT)
!      UMAX    = GLMAX(U3,3)
!C
!      RETURN
!      END
c------------------------------------------------------------------
!      SUBROUTINE CUMAXA(V1,V2,V3,UMAX)
!C
!      INCLUDE 'SIZE'
!      INCLUDE 'WZ'
!      INCLUDE 'GEOM'
!      INCLUDE 'INPUT'
!      INCLUDE 'ADAPT'
!C
!      DIMENSION V1(LX1,LY1,LZ1,1)
!     $        , V2(LX1,LY1,LZ1,1)
!     $        , V3(LX1,LY1,LZ1,1)
!      real  vxe(lx1u*ly1u*lz1u),vye(lx1u*ly1u*lz1u),vze(lx1u*ly1u*lz1u)
!      real  u(lx1u*ly1u*lz1u),v(lx1u*ly1u*lz1u),w(lx1u*ly1u*lz1u)
!      real  X(lx1u*ly1u*lz1u),R(lx1u*ly1u*lz1u)
!      DIMENSION U3(3)
!      integer nxu(ldimt1),nyu(ldimt1),nzu(ldimt1),ifld
!      INTEGER ICALLD
!      SAVE    ICALLD
!      DATA    ICALLD /0/
!C
!      NTOT  = NX1*NY1*NZ1
!      NTOTL = LX1U*LY1U*LZ1U
!      NTOTD = NTOTL*NDIM
!C
!C     Compute isoparametric partials.
!C
!     
!      DO IE=1,NELV
!        nxu = 0; nyu = 0; nzu = 0;
!        do ifld=1,ldimt
!          call getord(ie,j)
!          nxu(j) = nx1
!          nyu(j) = ny1
!          nzu(j) = nz1
!        end do
!        ifld = maxloc(nxu)
!        nx1 = nxu(ifld)
! 
!        CALL XYZRST1(XRM1E,YRM1E,ZRM1E,XSM1E,YSM1E,ZSM1E,XTM1E,YTM1E,
!     $             ZTM1E,IFAXIS)
!C   
!C     Compute maximum U/DX
!C
!!      IF (ICALLD.EQ.0) THEN
!!         ICALLD=1
!!         DRST (1)=ABS(ZGM1(2,1)-ZGM1(1,1))
!!         DRSTI(1)=1.0/DRST(1)
!!         DO 400 I=2,NX1-1
!!            DRST (I)=ABS(ZGM1(I+1,1)-ZGM1(I-1,1))/2.0
!!            DRSTI(I)=1.0/DRST(I)
!! 400     CONTINUE
!!         DRST (NX1)=DRST(1)
!!         DRSTI(NX1)=1.0/DRST(NX1)
!!      ENDIF
!C
!C     Zero out scratch arrays U,V,W for ALL declared elements...
!C
!        CALL RZERO3 (U,V,W,NTOTL)
!C
!        call MapV(vxe,vye,vze,ie,ifld)  
!        IF (NDIM.EQ.2) THEN
!
!        CALL VDOT2  (U,VXE,VYE,RXM1E,RYM1E,NTOT)
!        CALL VDOT2  (R,RXM1E,RYM1E,RXM1E,RYM1E,NTOT)
!        CALL VDOT2  (X,XRM1,YRM1,XRM1,YRM1,NTOT)
!        CALL COL2   (R,X,NTOT)
!        CALL VSQRT  (R,NTOT)
!        CALL INVCOL2(U,R,NTOT)
!C
!        CALL VDOT2  (V,VXE  ,VYE  ,SXM1E,SYM1E,NTOT)
!        CALL VDOT2  (R,SXM1E,SYM1E,SXM1E,SYM1E,NTOT)
!        CALL VDOT2  (X,XSM1E,YSM1E,XSM1E,YSM1E,NTOT)
!        CALL COL2   (R,X,NTOT)
!        CALL VSQRT  (R,NTOT)
!        CALL INVCOL2(V,R,NTOT)
!C
!        ELSE
!C
!        CALL VDOT3  (U,VXE,VYE,VZE,RXM1E,RYM1E,RZM1E,NTOT)
!        CALL VDOT3  (R,RXM1E,RYM1E,RZM1E,RXM1E,RYM1E,RZM1E,NTOT)
!        CALL VDOT3  (X,XRM1E,YRM1E,ZRM1E,XRM1E,YRM1E,ZRM1E,NTOT)
!        CALL COL2   (R,X,NTOT)
!        CALL VSQRT  (R,NTOT)
!        CALL INVCOL2(U,R,NTOT)
!C
!        CALL VDOT3  (V,VXE,VYE,VZE,SXM1E,SYM1E,SZM1E,NTOT)
!        CALL VDOT3  (R,SXM1E,SYM1E,SZM1E,SXM1E,SYM1E,SZM1E,NTOT)
!        CALL VDOT3  (X,XSM1E,YSM1E,ZSM1E,XSM1E,YSM1E,ZSM1E,NTOT)
!        CALL COL2   (R,X,NTOT)
!        CALL VSQRT  (R,NTOT)
!        CALL INVCOL2(V,R,NTOT)
!C
!        CALL VDOT3  (W,VXE,VYE,VZE,TXM1E,TYM1E,TZM1E,NTOT)
!        CALL VDOT3  (R,TXM1E,TYM1E,TZM1E,TXM1E,TYM1E,TZM1E,NTOT)
!        CALL VDOT3  (X,XTM1E,YTM1E,ZTM1E,XTM1E,YTM1E,ZTM1E,NTOT)
!        CALL COL2   (R,X,NTOT)
!        CALL VSQRT  (R,NTOT)
!        CALL INVCOL2(W,R,NTOT)
!C
!        ENDIF
!C
!        DO 500 IX=1,NX1
!        DO 500 IY=1,NY1
!        DO 500 IZ=1,NZ1
!              U()=ABS(U*DRSTI(IX) )
!              V()=ABS(V()*DRSTI(IY) )
!              W()=ABS(W()*DRSTI(IZ) )
!  500      CONTINUE
!      end do
!
!      U3(1)   = VLMAX(U,NTOT)
!      U3(2)   = VLMAX(V,NTOT)
!      U3(3)   = VLMAX(W,NTOT)
!      UMAX    = GLMAX(U3,3)
!C
!      RETURN
!      END
c------------------------------------------------------------------
      subroutine ComputeDxyz

      include "SIZE"
      include "DXYZ"
      include "WZ"
      include "ADAPT" 
      include "GEOM"
 
 
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      COMMON /FASTAXA/ WDDXA(LX1U*LY1U),WDDYTA(LY1U*LY1U),
     $                 WDDZTA(LZ1U*LZ1U)
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      real zgm1a(lx1u*3)
      integer ctr,ind1,ind2,ptr,ptr1,ptr2
 
      if (ndim.eq.2) then
        CALL ZWGLL (ZGM1A(1),WXM1E,NX1)
        CALL ZWGLL (ZGM1A(NX1+1),WYM1E,NY1)
        ZGM1(NZ1,3) = 0.
        WZM1(NZ1)   = 1.
        CALL DGLL (DXM1E,DXTM1E,ZGM1A(1),NX1,NX1)
        CALL DGLL (DYM1E,DYTM1E,ZGM1A(NX1+1),NY1,NY1)
        ptr = 0 
        DO 100 IY=1,NY1
        DO 100 IX=1,NX1
        ptr = ptr+1
        W3M1E(ptr)=WXM1E(IX)*WYM1E(IY)
  100   CONTINUE
      else
        CALL ZWGLL (ZGM1A(1),WXM1E,NX1)
        CALL ZWGLL (ZGM1A(NX1+1),WYM1E,NY1)
        CALL ZWGLL (ZGM1A(NX1+NY1+1),WZM1E,NZ1)
        CALL DGLL (DXM1E,DXTM1E,ZGM1A(1),NX1,NX1)
        CALL DGLL (DYM1E,DYTM1E,ZGM1A(NX1+1),NY1,NY1)
        CALL DGLL (DZM1E,DZTM1E,ZGM1A(NX1+NY1+1),NZ1,NZ1)
        ptr = 0
        DO 700 IZ=1,NZ1
        DO 700 IY=1,NY1
        DO 700 IX=1,NX1
        ptr = ptr+1
        W3M1E(ptr)=WXM1E(IX)*WYM1E(IY)*WZM1E(IZ)
  700   CONTINUE
      end if


      NXX = NX1*NX1
      CALL RZERO(WDDXA,NXX)
      DO 200 I=1,NX1
      DO 200 J=1,NX1
         ctr = (J-1)*NX1+I
      DO 200 IP=1,NX1
         ind1 = (I-1)*NX1+IP
         ind2 = (J-1)*NX1+IP
         WDDXA(ctr) = WDDXA(ctr) + WXM1E(IP)*DXM1E(ind1)*DXM1E(ind2)
  200 CONTINUE
      NYY=NY1*NY1
      CALL RZERO(WDDYTA,NYY)
      ctr = 0
      DO 300 I=1,NY1
      DO 300 J=1,NY1
         ctr = ctr+1
      DO 300 IP=1,NY1
         ind1 = (I-1)*NX1+IP
         ind2 = (J-1)*NX1+IP
         WDDYTA(ctr) = WDDYTA(ctr) + WYM1E(IP)*DYM1E(ind1)*DYM1E(ind2)
  300 CONTINUE
      NZZ=NZ1*NZ1
      CALL RZERO(WDDZTA,NZZ)
      ctr = 0
      DO 400 I=1,NZ1
      DO 400 J=1,NZ1
         ctr = ctr + 1
      DO 400 IP=1,NZ1
         ind1 = (I-1)*NX1+IP
         ind2 = (J-1)*NX1+IP
         WDDZTA(ctr) = WDDZTA(ctr) + WZM1E(IP)*DZM1E(ind1)*DZM1E(ind2)
  400 CONTINUE
 

!      IF (NDIM.EQ.3) THEN
!         DO 1001 IE=1,NELT
!            IF (.NOT.IFDFRM(IE)) THEN
!               DO 1000 IZ=1,NZ1
!               DO 1000 IY=1,NY1
!               DO 1000 IX=1,NX1
!                  G4M1(IX,IY,IZ,IE)=G1M1(IX,IY,IZ,IE)/WXM1E(IX)
!                  G5M1(IX,IY,IZ,IE)=G2M1(IX,IY,IZ,IE)/WYM1E(IY)
!                  G6M1(IX,IY,IZ,IE)=G3M1(IX,IY,IZ,IE)/WZM1E(IZ)
! 1000          CONTINUE
!            ENDIF
! 1001    CONTINUE
!      ELSE
!         DO 2001 IE=1,NELT
!            IF (.NOT.IFDFRM(IE)) THEN
!               DO 2000 IY=1,NY1
!               DO 2000 IX=1,NX1
!                  G4M1(IX,IY,1,IE)=G1M1(IX,IY,1,IE)/WXM1E(IX)
!                  G5M1(IX,IY,1,IE)=G2M1(IX,IY,1,IE)/WYM1E(IY)
! 2000          CONTINUE
!            ENDIF
! 2001    CONTINUE
!      ENDIF

      return
      end 
c------------------------------------------------------------------
      subroutine UpdateSol

      include 'SIZE'
      include 'ADAPT'
      include 'INPUT'

      integer ifld,ntot,ptr,ptr1,ptr2

      !if (ifflow) then
      !  ifld = 1
      !  ntot = ntota(ifld)
      !  call copy(soltmp,vxa,ntot)
      !  call interp(vxa,soltmp,ifld)
      !  call copy(soltmp,vya,ntot)
      !  call interp(vya,soltmp,ifld)
      !  if (ndim.eq.3) then
      !    call copy(soltmp,vza,ntot)
      !    call interp(vza,soltmp,ifld)
      !  end if
      !  !call copy(soltmp,vxlaga,ntot)
      !  !do i=1,2   ! Always 2 lag velocities ???
      !  !  ptr1 = 
      !  !  ptr2 = 
      !  !  call interp(vxlaga(ptr2),soltmp,ifld)
      !  !end do
      !  !call copy(soltmp,vylaga,ntot)
      !  !do i=1,2   
      !  !  ptr1 = 
      !  !  ptr2 = 
      !  !  call interp(vylaga(ptr2),soltmp(ptr1),ifld)
      !  !end do
      !  !if (ndim.eq.3) then
      !  !  call copy(soltmp,vzlaga,ntot)
      !  !  do i=1,2   
      !  !    ptr1 = 
      !  !    ptr2 = 
      !  !    call interp(vzlaga(ptr2),soltmp(ptr1),ifld)
      !  !  end do
      !  !end if
      !end if

       ! Scalar solutions
         if (nid==1) then
       !print *, 'Update Sol'
       !print '(16F8.3)', tad
       !print *
        end if
       ntot = ivlsum(ntota(2),ldimt)
       call copy(soltmp,tad,ntot)
       do ifld=2,ldimt1
           ptr1 = adptr(1,ifld)
           ptr2 = newptr(1,ifld)
           call interp(tad(ptr2),soltmp(ptr1),ifld)
       end do
         if (nid==1) then
       !print '(25F8.3)', tad
       !print *  
         end if 
       ! Scalar lag solutions
       !print *, 'Update Lag Start'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tlaga
       !elseif (nx1.eq.4) then
       !  print '(25F8.3)', tlaga
       !else
       !  print '(25F8.3)', tlaga
       !end if
       !print *
       ntot = ntot*(lorder-1)
       call copy(soltmp,tlaga,ntot)
      do ifld=2,ldimt1
       if (adaptflag(1,ifld)) then
         do i=1,lorder-1
           ptr1=ivlsum(ntota(2),(ifld-2))*(lorder-1)+ntota(ifld)*(i-1)+1
           ptr2=ivlsum(newntota(2),(ifld-2))*(lorder-1)
     $           + newntota(ifld)*(i-1) + 1
          ! print *, 'UPDATE LAG PTRS', ptr1, ptr2
           call interp(tlaga(ptr2),soltmp(ptr1),ifld)
          ! print *, ifld, i, ptr1, ptr2
          ! if (nx1.eq.3) then
          !   print *, 'soltmp'
          !   print '(25F8.3)', soltmp
          !   print *, 'tlaga'
          !   print '(25F8.3)', tlaga
          ! elseif (nx1.eq.4) then
          !   print *, 'soltmp'
          !   print '(25F8.3)', soltmp
          !   print *, 'tlaga'
          !   print '(25F8.3)', tlaga
          ! else
          !   print *, 'soltmp'
          !   print '(25F8.3)', soltmp
          !   print *, 'tlaga'
          !   print '(25F8.3)', tlaga
          ! end if
         end do
       end if
      end do

      ! Update vgradt1a and vgradt2a
      ntot = ivlsum(ntota(2),ldimt)
      call copy(soltmp,vgradt1a,ntot)
      ptr1 = adptr(1,2)
      call col2(soltmp,binvm1a(ptr1+ntota(1)),ntot)

      !print *, 'VGRADT1A'
      !print '(25F8.3)', vgradt1a
      !print *

      do ifld=2,ldimt1
          ptr1 = adptr(1,ifld)
          ptr2 = newptr(1,ifld)
          call interp(vgradt1a(ptr2),soltmp(ptr1),ifld)
      end do
      !print *, 'VGRADT1A'
      !print '(25F8.3)', vgradt1a
      !print *
      call copy(soltmp,vgradt2a,ntot)
      ptr1 = adptr(1,2)
      call col2(soltmp,binvm1a(ptr1+ntota(1)),ntot)
      !print *, 'VGRADT2A'
      !print '(25F8.3)', vgradt2a
      !print *
      do ifld=2,ldimt1
          ptr1 = adptr(1,ifld)
          ptr2 = newptr(1,ifld)
          call interp(vgradt2a(ptr2),soltmp(ptr1),ifld)
      end do
      !print *, 'VGRADT2A'
      !print '(25F8.3)', vgradt2a
      !print *

       !print *, 'Update Lag End'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tlaga
       !else
       !  print '(25F8.3)', tlaga
       !end if
       !print *

      end
c------------------------------------------------------------------
      subroutine UpdateMass

      ! Update mass matrices for each field     

      include 'SIZE'
      include 'ADAPT'
      include 'WZ'
      include 'INPUT'
      include 'TSTEP'

      real wk1(lx1u),wk2(lx1u*lx1u),wka(lx1u*ly1u*lz1u),
     $     wkb(lx1u*ly1u*lz1u),wkc(lx1u*ly1u*lz1u),wk3(lx1u*ly1u*lz1u),
     $     weight
      integer iord,ctr,ntot,i,nxyz,nxyz1,ifld,bptr1,bptr2,ptr,ptr1,ptr2,
     $        bptr3,bptr4
      real  XRM1E(LX1U*LY1U*LZ1U),YRM1E(LX1U*LY1U*LZ1U),
     $      XSM1E(LX1U*LY1U*LZ1U),YSM1E(LX1U*LY1U*LZ1U),
     $      XTM1E(LX1U*LY1U*LZ1U),YTM1E(LX1U*LY1U*LZ1U),
     $      ZRM1E(LX1U*LY1U*LZ1U),ZSM1E(LX1U*LY1U*LZ1U),
     $      ZTM1E(LX1U*LY1U*LZ1U)

      ntot = ivlsum(ntota,ldimt1)
      call copy(soltmp,bm1a,ntot)
      call copy(soltmpb,areaa,lx1u*lz1u*lelt*ldimt1*2*ldim) 
      do ifld=1,ldimt1
        ntot = ntota(ifld)
        nel = nelfld(ifld) 
        do ie=1,nel
         ptr1 = adptr(ie,ifld) 
         ptr2 = newptr(ie,ifld)
         if (ifld.ne.1) then
            ptr1 = ptr1+ntota(1)
            ptr2 = ptr2+newntota(1)
         endif
         bptr1 = bptr(ie,ifld)
         bptr2 = newbptr(ie,ifld)
         if (adaptflag(ie,ifld)) then
           call getord(ie,ifld) 
           nx1 = nx1+adinc; ny1 = ny1+adinc;
           if (nz1.ne.1) nz1 = nz1+adinc
           nxyz = nx1*ny1*nz1
           ! Get weights
           CALL ZWGLL (wk1,WXM1E,NX1)
           CALL ZWGLL (wk1,WYM1E,NY1)
           if (ndim.eq.2) then
              do i=1,nz1
                WZM1E(i) = 1.
              end do
           else
             CALL ZWGLL (wk1,WZM1E,NZ1)
           end if
           ctr = 1
           DO 100 IZ=1,NZ1
           DO 100 IY=1,NY1
           DO 100 IX=1,NX1
           W3M1E(ctr)=WXM1E(IX)*WYM1E(IY)*WZM1E(IZ)
           ctr = ctr+1
  100      CONTINUE
           !! Get Jacobian
           ! Get x, y, z for this element
           call genxyz1(xm1e,ym1e,zm1e,nx1,ny1,nz1,ie)
           
           NXY1  = NX1*NY1
           NYZ1  = NY1*NZ1
           NXZ1  = NX1*NZ1
           NXYZ1 = NX1*NY1*NZ1
           ntot1 = nxyz1
           ! Get new derivative matrices
           call ComputeDxyz 
      
           ! Get partial derivatives for this element
           CALL XYZRST1(XRM1E,YRM1E,ZRM1E,XSM1E,YSM1E,ZSM1E,XTM1E,
     $                  YTM1E,ZTM1E,IFAXIS,IE)
                      
           IF (NDIM.EQ.2) THEN
             CALL RZERO   (JACM1E,NTOT1)
             CALL ADDCOL3 (JACM1E,XRM1E,YSM1E,NTOT1)
             CALL SUBCOL3 (JACM1E,XSM1E,YRM1E,NTOT1)
              

           ELSE
             CALL RZERO   (JACM1E,NTOT1)
             CALL ADDCOL4 (JACM1E,XRM1E,YSM1E,ZTM1E,NTOT1)
             CALL ADDCOL4 (JACM1E,XTM1E,YRM1E,ZSM1E,NTOT1)
             CALL ADDCOL4 (JACM1E,XSM1E,YTM1E,ZRM1E,NTOT1)
             CALL SUBCOL4 (JACM1E,XRM1E,YTM1E,ZSM1E,NTOT1)
             CALL SUBCOL4 (JACM1E,XSM1E,YRM1E,ZTM1E,NTOT1)
             CALL SUBCOL4 (JACM1E,XTM1E,YSM1E,ZRM1E,NTOT1)

             ! AREAA
             
             CALL VCROSS(wka,wkb,wkc,XSM1E,YSM1E,ZSM1E,XTM1E,YTM1E,
     $                   ZTM1E,NTOT1)
             CALL VDOT3 (wk3,wka,wkb,wkc,wka,wkb,wkc,NTOT1)
             !print *, 'NX1'
             !print *, nx1,ny1,nz1,ntot1
             do 200 iz=1,nz1
               do 200 iy=1,ny1
                 weight = wym1e(iy)*wzm1e(iz)
                 bptr3 = bptr2+iy+(iz-1)*ny1+nyz1-1
                 bptr4 = bptr2+iy+(iz-1)*ny1+nyz1*3-1
                 ind1 = (iz-1)*nxy1+(iy-1)*nx1+nx1
                 ind2 = (iz-1)*nxy1+(iy-1)*nx1+1
                 areaa(bptr3) = sqrt(wk3(ind1))*weight
                 areaa(bptr4) = sqrt(wk3(ind2))*weight
!               print *, 'BPTRS 1',bptr1,bptr3,bptr4,
!     $          bptr3-bptr1+1, bptr4-bptr1+1,areaa(bptr3),areaa(bptr4)
  200        continue
             
             CALL VCROSS(wka,wkb,wkc,XRM1E,YRM1E,ZRM1E,XTM1E,YTM1E,
     $                   ZTM1E,NTOT1)
             CALL VDOT3 (wk3,wka,wkb,wkc,wka,wkb,wkc,NTOT1)
             do 300 iz=1,nz1
               do 300 ix=1,nx1
                 weight = wxm1e(ix)*wzm1e(iz)
                 bptr3 = bptr2+ix+(iz-1)*nx1-1
                 bptr4 = bptr2+ix+(iz-1)*nx1+nxz1*2-1
                 ind1 = (iz-1)*nxy1+ix
                 ind2 = (iz-1)*nxy1+(ny1-1)*nx1+ix
                 areaa(bptr3) = sqrt(wk3(ind1))*weight
                 areaa(bptr4) = sqrt(wk3(ind2))*weight
!               print *, 'BPTRS 2',bptr1,bptr3,bptr4,
!     $          bptr3-bptr1+1, bptr4-bptr1+1,areaa(bptr3),areaa(bptr4)
  300        continue
             
             CALL VCROSS(wka,wkb,wkc,XRM1E,YRM1E,ZRM1E,XSM1E,YSM1E,
     $                   ZSM1E,NTOT1)
             CALL VDOT3 (wk3,wka,wkb,wkc,wka,wkb,wkc,NTOT1)
             do 400 ix=1,nx1
               do 400 iy=1,ny1
                 weight = wxm1e(ix)*wym1e(iy)
                 bptr3 = bptr2+iy+(ix-1)*ny1+nxy1*4-1
                 bptr4 = bptr2+iy+(ix-1)*ny1+nxy1*5-1
                 ind1 = (iy-1)*nx1+ix
                 ind2 = (nz1-1)*nxy1+(iy-1)*nx1+ix
                 areaa(bptr3) = sqrt(wk3(ind1))*weight
                 areaa(bptr4) = sqrt(wk3(ind2))*weight
!               print *, 'BPTRS 3',bptr1,bptr3,bptr4,
!     $          bptr3-bptr1+1, bptr4-bptr1+1,areaa(bptr3),areaa(bptr4)
  400        continue

           END IF
           ! Build elemental mass matrix
           nxyz = nx1*ny1*nz1
           CALL COL3 (BM1A(ptr2),JACM1E,W3M1E,NXYZ)
         else
           nxyz = nx1*ny1*nz1
           call copy(bm1a(ptr2),soltmp(ptr1),nxyz)
           nxyz2 = nx1*nz1*2*ldim
           call copy(areaa(ptr2),soltmpb(ptr1),nxyz2)
         end if
        end do
      end do

      ! Update inverse mass matrix
      do i=1,ldimt1
        ntot = newntota(i)
        ptr = newptr(1,i)
        if (i.ne.1) ptr = ptr + newntota(1)
        call getord(1,i)
        if (adaptflag(1,i)) then
          nx1 = nx1 + adinc
          ny1 = ny1 + adinc
          if (nz1.ne.1) nz1 = nz1+adinc
        end if
        call copy(binvm1a(ptr),bm1a(ptr),ntot)
        call dssum(binvm1a(ptr),nx1,ny1,nz1)
        call invcol1(binvm1a(ptr),ntot) 
      end do 

      return
      end
c------------------------------------------------------------------
c      subroutine UpdateMask
c
c
c      ! scalar masks
c      ntot = ivlsum(ntota(2),ntot,ldimt) 
c      call copy(soltmp,tmask,ntot)
c      do 80 ifld=2,ldimt1
c        do 80 ie=1,nelt
c          ptr1 = adptr(ie,ifld)
c          ptr2 = newptr(ie,ifld)
c          call getord(ie,ifld) 
c          if (adaptflag(ie,ifld)) then
c            nx1 = nx1+1; ny1 = ny1+1;
c            if (nz1.ne.1) nz1 = nz1+1
c            ! Set interior points to one
c            nxy = nx1*ny1
c            do 90 k=2,nz1-1
c              do 90 j=2,ny1-1
c                do 90 i=2,nx1-1
c                  ind = ptr2+(k-1)*nxy+(j-1)*nx1+i-1
c                  tmask(ind) = 1.0  
c 90         continue
c            ! 
c          else
c            nxyz = nx1*ny1*nz1
c            call copy(tmask(ptr2),soltmp(ptr1),nxyz)
c          end if
c 80   continue
c
c      end 
c------------------------------------------------------------------
      subroutine UpdateMult
 
      include 'SIZE'
      include 'ADAPT'
      include 'TSTEP'
      include 'INPUT'
      
      integer ntot,ptr,ptr1,ptr2

      if (ifflow) then
         ifield = 1
         ntot = ntota (ifield)
         call getord  (1,ifield) ! Uniform order for now
         call rone    (vmulta,ntot)
         call dssum   (vmulta,nx1,ny1,nz1)
         call invcol1 (vmulta,ntot)
      endif
      if (ifheat) then
         ifield = 2
         ntot = ntota (ifield)
         call getord  (1,ifield)
         call rone    (tmulta,ntot)
         call dssum(tmulta,nx1,ny1,nz1)
         call invcol1 (tmulta,ntot)
      endif
      if (.not.ifflow) call copy(vmulta,tmulta,ntot)
      do ifield=3,nfield                  ! Additional pass. scalrs.
         ntot = ntota (ifield)
         ptr1 = adptr (1,ifield)
         call getord  (1,ifield)
         call rone    (tmulta(ptr1),ntot)
         call dssum(tmulta(ptr1),nx1,ny1,nz1)
         call invcol1 (tmulta(ptr1),ntot)
      enddo
    
      return
      end 
c------------------------------------------------------------------
      subroutine UpdateBCs
        
      include 'SIZE'
      include 'ADAPT'
      include 'INPUT'
      include 'TSTEP'

      common /nekcb/ cb
      integer ifld,iel,i,ntot,iface,nfaces,ptr,ptr1,ptr2
      character*3 cb
      character*1 cb1(3)
      equivalence (cb1,cb)

      ! Update boundary masks
      !call bcmask
      NFACES=2*NDIM
 
      IF (IFHEAT) THEN
        ntot = ntota(ifield)
        DO 1200 IFIELD=2,NFIELD
           ptr1 = adptr(1,ifield)
           NEL    = NELFLD(IFIELD)
           NTOT   = ntota(ifield)
           CALL RONE (TMASKA(ptr1),NTOT)
        DO 1100 IEL=1,NEL
           call getord(iel,ifield)
           ptr2 = adptr(iel,ifield)
        DO 1100 IFACE=1,NFACES
           CB =CBC(IFACE,IEL,IFIELD)
C
C          Assign mask values.
C
           IF  (CB.EQ.'T  ' .OR. CB.EQ.'t  ' .OR. 
     $          CB.EQ.'MCI' .OR. CB.EQ.'MLI' .OR.
     $          CB.EQ.'KD ' .OR. CB.EQ.'kd ' .OR.
     $          CB.EQ.'ED ' .OR. CB.EQ.'ed ' .OR.
     $          CB.EQ.'KW ' .OR. CB.EQ.'KWS' .OR. CB.EQ.'EWS')
     $          CALL FACEV (TMASKA(ptr2),
     $                      1,IFACE,0.0,NX1,NY1,NZ1)
 1100      CONTINUE
        CALL DSOP (TMASKA(ptr1),'MUL',NX1,NY1,NZ1)
 1200   CONTINUE
      ENDIF
 


      ! Update boundaries for adapted solutions (including lags)
      if (ifflow) then
        ! Update flow boundaries

      end if

       !print *, 'Update T'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tad
       !elseif (nx1.eq.4) then
       !  print '(25F8.3)', tad
       !else
       !  print '(25F8.3)', tad
       !end if
       !print *

       !print *, 'Update Lag'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tlaga
       !elseif (nx1.eq.4) then
       !  print '(25F8.3)', tlaga
       !else
       !  print '(25F8.3)', tlaga
       !end if
       !print *
 
      ! Update scalars
      do ifld=2,ldimt1
        if (adaptflag(1,ifld)) then
          ifield = ifld
          ptr = adptr(1,ifld)   ! Uniform order 
          call BCDIRSCA(Tad(ptr))
          do i=1,lorder-1
          ptr = ivlsum(ntota(2),(ifld-2))*(lorder-1)+ntota(ifld)*(i-1)+1
            call BCDIRSCA(tlaga(ptr))
          end do
        end if
      end do
       !print *, 'Update T'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tad
       !elseif (nx1.eq.4) then
       !  print '(25F8.3)', tad
       !else
       !  print '(25F8.3)', tad
       !end if
       !print *
       !print *, 'Update Lag'
       !if (nx1.eq.3) then
       !  print '(25F8.3)', tlaga
       !elseif (nx1.eq.4) then
       !  print '(25F8.3)', tlaga
       !else
       !  print '(25F8.3)', tlaga
       !end if
       !print *

      return
      end
c------------------------------------------------------------------
      subroutine cdscala(igeom)
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      LOGICAL          IFCONV
C
      COMMON /SCRNSA/ TA(LX1U*LY1U*LZ1U*LELT)
     $              ,TB(LX1U*LY1U*LZ1U*LELT)
      COMMON /SCRVHA/ H1(LX1U*LY1U*LZ1U*LELT)
     $              ,H2(LX1U*LY1U*LZ1U*LELT)

      integer ptra,ptr,ptr1,ptr2

      include 'ORTHOT'
      include 'ADAPT'

      napprox(1) = laxt

      nel    = nelfld(ifield)
      ntot = ntota(ifield) 

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         ! Need to change these to account for variable order
         if (nid==0) then
         print *, 'BQ', igeom, 1, nid
         print '(25F8.3)', BQA
         print *
         end if
         call makeqa
         if (nid==0) then
         print *, 'BQ', igeom, 2, nid
         print '(25F8.3)', BQA
         print *
         end if
         call lagscal

      else                   ! geometry at t^n

         IF (IFPRINT) THEN
         IF (IFMODEL .AND. IFKEPS) THEN
            NFLDT = NFIELD - 1
            IF (IFIELD.EQ.NFLDT.AND.NID.EQ.0) THEN
               WRITE (6,*) ' Turbulence Model - k/epsilon solution'
            ENDIF
         ELSE
            IF (IFIELD.EQ.2.AND.NID.EQ.0) THEN
               WRITE (6,*) ' Temperature/Passive scalar solution'
            ENDIF
         ENDIF
         ENDIF
         if1=ifield-1
         write(name4,1) if1-1
    1    format('PS',i2)
         if(ifield.eq.2) write(name4,'(A4)') 'TEMP'

C
C        New geometry
C
         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)
         INTYPE = 0
         IF (IFTRAN) INTYPE = -1
         ptra = adptr(1,ifield)
         CALL SETHLMA (H1,H2,INTYPE)   
         CALL BCNEUSCA(TA,-1)
         CALL ADD2    (H2,TA,NTOT) 
         CALL BCDIRSCA(TAD(ptra))
         !if (nid==0) then
         !print *, 'AX PARAMS', nid 
         !print *, 'T', nid
         !print '(16F8.3)', TAD
         !print *
         !print *, 'H1', nid
         !print '(16F8.3)', H1
         !print *
         !print *, 'H2', nid
         !print '(16F8.3)', H2
         !print *
         !end if
         CALL AXHELMA (TA,TAD(ptra),H1,H2,IMESH,isd)
         CALL SUB3    (TB,BQA(ptra),TA,NTOT)
         !if (nid==0) then
         !print *, 'first', nid
         !print *, 'BQ', nid
         !print '(16F8.3)', BQA
         !print *
         !print *, 'TA', nid
         !print '(16F8.3)', TA
         !print *
         !print *, 'TB', nid
         !print '(16F8.3)', TB
         !print *
         !end if
         CALL BCNEUSCA(TA,1)
         CALL ADD2    (TB,TA,NTOT)

         !if (nid==0) then
         !print *, 'TA',nid
         !print '(16F8.3)', TA
         !print *
         !print *, 'TB',nid
         !print '(16F8.3)', TB
         !print *
         !print *, 'H1',nid
         !print '(16F8.3)', H1
         !print *
         !print *, 'H2',nid
         !print '(16F8.3)', H2
         !print *
         !print *, 'tmask', ptra, nid
         !print '(16F8.3)', tmaska
         !print *
         !print *, 'tmult', nid
         !print '(16F8.3)', tmulta
         !print *
         !print *, 'imesh, tolht, nmxh', nid
         !print '(I8,F8.3,I8)', imesh, tolht(ifield),nmxh
         !print *
         !end if
        CALL HMHOLTZA(name4,TA,TB,H1,H2
     $                 ,TMASKA(ptra),TMULTA(ptra)
     $                 ,IMESH,TOLHT(IFIELD),NMXH,1)

         !if(iftmsh(ifield)) then
         !  call hsolve  (name4,TA,TB,H1,H2 
     $   !                ,tmask(ptr)
     $   !                ,tmult(ptr)
     $   !                ,imesh,tolht(ifield),nmxh,1
     $   !                ,approx,napprox,bintm1)
         !else
         !  call hsolve  (name4,TA,TB,H1,H2 
     $   !                ,tmask(ptr)
     $   !                ,tmult(ptr)
     $   !                ,imesh,tolht(ifield),nmxh,1
     $   !                ,approx,napprox,binvm1)
         !endif 

        ! if (nid==0) then
         !print *,'TA',nid
         !print '(16F8.3)', ta
         !print *
        ! end if
         call add2    (tad(ptra),ta,ntot)
        ! if (nid==0) then
         !print *,'T',nid
         !print '(16F8.3)', tad
         !print *
        ! end if
         call cvgnlps (ifconv) ! Check convergence for nonlinear problem 
         if (ifconv) goto 2000

C        Radiation case, smooth convergence, avoid flip-flop (ER).
         CALL CMULT (TA,0.5,NTOT)
         CALL SUB2  (TAD(ptra),TA,NTOT)

 1000    CONTINUE
 2000    CONTINUE
         CALL BCNEUSCA(TA,1)
         CALL ADD2 (BQA(ptra),TA,NTOT) ! no idea why... pf
      endif

 
      return
      end
c=======================================================================
      subroutine axhelma(au,u,helm1,helm2,imesh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ADAPT'
C
      COMMON /FASTAXA/ WDDXA(LX1U*LX1U),WDDYTA(LY1U*LY1U),
     $                 WDDZTA(LZ1U*LZ1U)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      COMMON /ISTEP2/ IFIELD
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV, GEOFLAG
C
      REAL           AU    (LX1U*LY1U*LZ1U)
     $ ,             U     (LX1U*LY1U*LZ1U)
     $ ,             HELM1 (LX1U*LY1U*LZ1U)
     $ ,             HELM2 (LX1U*LY1U*LZ1U)
      COMMON /CTMP1/ DUDR  (LX1U*LY1U*LZ1U)
     $ ,             DUDS  (LX1U*LY1U*LZ1U)
     $ ,             DUDT  (LX1U*LY1U*LZ1U)
     $ ,             TMP1  (LX1U*LY1U*LZ1U)
     $ ,             TMP2  (LX1U*LY1U*LZ1U)
     $ ,             TMP3  (LX1U*LY1U*LZ1U)

      REAL           TM1   (LX1U*LY1U*LZ1U)
      REAL           TM2   (LX1U*LY1U*LZ1U)
      REAL           TM3   (LX1U*LY1U*LZ1U)
      REAL           DUAX  (LX1U)
      REAL           YSM1  (LX1U)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)
      INTEGER        ptr,ptr1,ptr2

C
      IF(IMESH.EQ.1) NEL=NELV
      IF(IMESH.EQ.2) NEL=NELT
      NXY=NX1*NY1
      NYZ=NY1*NZ1
      NXZ=NX1*NZ1
      NXYZ=NX1*NY1*NZ1
      ntot = ntota(ifield)  ! KED 
C
#ifndef NOTIMER
      if (icalld.eq.0) taxhm=0.0
      icalld=icalld+1
      naxhm=icalld
      etime1=dnekclock()
#endif
C
      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
      CALL RZERO (AU,NTOT)
C
      !  print *, 'AX VALS'
      DO 100 IEL=1,NEL
        ie=iel
        nx1o = nx1
        if (iel.eq.1) nx1o = 0
        call getord(ie,ifield)
        NXY=NX1*NY1
        NYZ=NY1*NZ1
        NXZ=NX1*NZ1
        NXYZ=NX1*NY1*NZ1
        geoflag = .false.
        !if (nx1.ne.lx1) then
          geoflag = .true. 
          call ComputeG(ie,GM1A)
        !else
        !  if (nx1o.eq.nx1) call ComputeDxyz 
        !end if

C
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
C
        IF (NDIM.EQ.2) THEN
C
C       2-d case ...............
C
           IF (IFFAST(IEL)) THEN
C
C          Fast 2-d mode: constant properties and undeformed element
C
           ptr1 = adptr1(ie,ifield)
           H1 = HELM1(ptr1)
           !print *, 'WDDX'
           !print '(16F8.3)', WDDXA
           !print *
           !print *, 'U', ptr1
           !print '(16F8.3)', U
           !print *
           !print *, 'nx1,ny1,nz1,nyz',nx1,ny1,nz1,nyz
           !print *
           !print *, 'WDDYT'
           !print '(16F8.3)', WDDYTA
           !print *
           CALL MXM   (WDDXA,NX1,U(ptr1),NX1,TM1,NYZ) ! Need to alter WDDX
           CALL MXM   (U(ptr1),NX1,WDDYTA,NY1,TM2,NY1)
           !print *, 'TM1'
           !print '(16F8.3)', TM1
           !print *
           !print *, 'TM2'
           !print '(16F8.3)', TM2
           !print *
           if (geoflag) then
             CALL COL2  (TM1,GM1A(nxyz*3+1),NXYZ)   
             CALL COL2  (TM2,GM1A(nxyz*4+1),NXYZ)  
           else
             CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ) 
             CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)  
           end if
           !print *, 'G', nxyz*3+1, nxyz*4+1, geoflag
           !print '(16F8.3)', GM1A
           !print *
           !print *, 'TM1 AGAIN'
           !print '(16F8.3)', TM1
           !print *
           !print *, 'TM2 AGAIN'
           !print '(16F8.3)', TM2
           !print *
           CALL ADD3  (AU(ptr1),TM1,TM2,NXYZ)
           !print *, 'AU', ptr1
           !print '(16F8.3)', AU
           !print *
           CALL CMULT (AU(ptr1),H1,NXYZ)
           !print *, 'AU again', ptr1
           !print '(16F8.3)', AU
           !print *
          
           ELSE
C
C          General case, speed-up for undeformed elements
C
           ptr1 = adptr1(ie,ifield)
           CALL MXM  (DXM1E,NX1,U(ptr1),NX1,DUDR,NYZ)
           CALL MXM  (U(ptr1),NX1,DYTM1E,NY1,DUDS,NY1)
           if (geoflag) then
             CALL COL3 (TMP1,DUDR,GM1A(1),NXYZ)
             CALL COL3 (TMP2,DUDS,GM1A(nxyz+1),NXYZ)
           else
             CALL COL3 (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
             CALL COL3 (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
           end if
           IF (IFDFRM(IEL)) THEN
             if (geoflag) then
               CALL ADDCOL3 (TMP1,DUDS,GM1A(3*nxyz+1),NXYZ)
               CALL ADDCOL3 (TMP2,DUDR,GM1A(3*nxyz+1),NXYZ)
             else              
               CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
               CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
             end if
           ENDIF
           CALL COL2 (TMP1,HELM1(ptr1),NXYZ)
           CALL COL2 (TMP2,HELM1(ptr1),NXYZ)
           CALL MXM  (DXTM1E,NX1,TMP1,NX1,TM1,NYZ)
           CALL MXM  (TMP2,NX1,DYM1E,NY1,TM2,NY1)
           CALL ADD2 (AU(ptr1),TM1,NXYZ)
           CALL ADD2 (AU(ptr1),TM2,NXYZ)
C
        ENDIF
C
        ELSE
C
C       3-d case ...............
C
           IF (IFFAST(IEL)) THEN
C
C          Fast 3-d mode: constant properties and undeformed element
C
           ptr1 = adptr1(ie,ifield)
           H1 = HELM1(ptr1)
           CALL MXM   (WDDXA,NX1,U(ptr1),NX1,TM1,NYZ)
           !print *, 'WDDXA'
           !print '(16F8.3)',wddxa
           !print *
           !print *, 'U', ptr1
           !print '(16F8.3)',U
           !print *
           !print *, 'nx1 nyz', nx1, nyz 
           !print *
           !print *, 'TM1'
           !print '(16F8.3)',tm1
           !print *
           DO 5 IZ=1,NZ1
           ptr2 = nxy*(IZ-1)
           CALL MXM   (U(ptr1+ptr2),NX1,WDDYTA,NY1,TM2(ptr2+1),NY1)
 5         CONTINUE
           !print *, 'TM2'
           !print '(16F8.3)',tm2
           !print *
           CALL MXM   (U(ptr1),NXY,WDDZTA,NZ1,TM3,NZ1)
           !print *, 'TM3'
           !print '(16F8.3)',tm3
           !print *
           if (geoflag) then
             CALL COL2  (TM1,GM1A(nxyz*3+1),NXYZ)
             CALL COL2  (TM2,GM1A(nxyz*4+1),NXYZ)
             CALL COL2  (TM3,GM1A(nxyz*5+1),NXYZ)
           else
             CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
             CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
             CALL COL2  (TM3,G6M1(1,1,1,IEL),NXYZ)
           end if
           !print *, 'TM1 AGAIN'
           !print '(16F8.3)',tm1
           !print *
           !print *, 'TM2 AGAIN'
           !print '(16F8.3)',tm2
           !print *
           !print *, 'TM3 AGAIN'
           !print '(16F8.3)',tm3
           !print *
           CALL ADD3  (AU(ptr1),TM1,TM2,NXYZ)
           CALL ADD2  (AU(ptr1),TM3,NXYZ)
           CALL CMULT (AU(ptr1),H1,NXYZ)
           !print *, 'AU'
           !print '(16F8.3)',au
           !print *
C
           ELSE
C
C          General case, speed-up for undeformed elements
C
           ptr1 = adptr1(ie,ifield)
           CALL MXM(DXM1E,NX1,U(ptr1),NX1,DUDR,NYZ)
           DO 10 IZ=1,NZ1
              ptr2 = (IZ-1)*nxy
              CALL MXM(U(ptr1+ptr2),NX1,DYTM1E,NY1,DUDS(ptr2+1),NY1)
   10      CONTINUE
           CALL MXM     (U(ptr1),NXY,DZTM1E,NZ1,DUDT,NZ1)
           if (geoflag) then
             CALL COL3    (TMP1,DUDR,GM1A(1),NXYZ)
             CALL COL3    (TMP2,DUDS,GM1A(nxyz+1),NXYZ)
             CALL COL3    (TMP3,DUDT,GM1A(2*nxyz+1),NXYZ)
           else
             CALL COL3    (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
             CALL COL3    (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
             CALL COL3    (TMP3,DUDT,G3M1(1,1,1,IEL),NXYZ)
           end if
           IF (IFDFRM(IEL)) THEN
              if (geoflag) then
                CALL ADDCOL3 (TMP1,DUDS,GM1A(3*nxyz+1),NXYZ)
                CALL ADDCOL3 (TMP1,DUDT,GM1A(4*nxyz+1),NXYZ)
                CALL ADDCOL3 (TMP2,DUDR,GM1A(3*nxyz+1),NXYZ)
                CALL ADDCOL3 (TMP2,DUDT,GM1A(5*nxyz+1),NXYZ)
                CALL ADDCOL3 (TMP3,DUDR,GM1A(4*nxyz+1),NXYZ)
                CALL ADDCOL3 (TMP3,DUDS,GM1A(5*nxyz+1),NXYZ)
              else
                CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
                CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
                CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
                CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
                CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
                CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
              end if
           ENDIF
           CALL COL2 (TMP1,HELM1(ptr1),NXYZ)
           CALL COL2 (TMP2,HELM1(ptr1),NXYZ)
           CALL COL2 (TMP3,HELM1(ptr1),NXYZ)
           CALL MXM  (DXTM1E,NX1,TMP1,NX1,TM1,NYZ)
           DO 20 IZ=1,NZ1
              ptr2 = (IZ-1)*nxy+1
              CALL MXM(TMP2(ptr2),NX1,DYM1E,NY1,TM2(ptr2),NY1)
   20      CONTINUE
           CALL MXM  (TMP3,NXY,DZM1E,NZ1,TM3,NZ1)
           CALL ADD2 (AU(ptr1),TM1,NXYZ)
           CALL ADD2 (AU(ptr1),TM2,NXYZ)
           CALL ADD2 (AU(ptr1),TM3,NXYZ)
C
           ENDIF
C
        ENDIF
C
 100  CONTINUE
     
      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr+ntota(1)

      ntot = ntota(ifield)

      IF (IFH2) CALL ADDCOL4 (AU,HELM2,BM1A(ptr),U,NTOT)
           !print *, 'AU FINAL',ntot
           !print '(16F8.3)', AU
           !print *
           !print *, 'BM', ptr
           !print '(16F8.3)', BM1A
           !print *
           !print *, 'U'
           !print '(16F8.3)', U
           !print *
           !print *, 'HELM2'
           !print '(16F8.3)', HELM2
           !print *

C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
!      IF (IFAXIS.AND.(ISD.EQ.2)) THEN
!         DO 200 IEL=1,NEL
!C
!            IF (IFRZER(IEL)) THEN
!               CALL MXM(U  (1,1,1,IEL),NX1,DATM1,NY1,DUAX,1)
!               CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
!            ENDIF
!C
!            DO 190 J=1,NY1
!            DO 190 I=1,NX1
!               IF (YM1(I,J,1,IEL).NE.0.) THEN
!                  TERM1 = BM1(I,J,1,IEL)*U(I,J,1,IEL)/YM1(I,J,1,IEL)**2
!                  IF (IFRZER(IEL)) THEN
!                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)*DUAX(I)
!     $                       *JACM1(I,1,1,IEL)/YSM1(I)
!                  ELSE
!                     TERM2 = 0.
!                  ENDIF
!                  AU(I,J,1,IEL) = AU(I,J,1,IEL)
!     $                          + HELM1(I,J,1,IEL)*(TERM1+TERM2)
!               ENDIF
!  190       CONTINUE
!  200    CONTINUE
!      ENDIF
C
#ifndef NOTIMER
      taxhm=taxhm+(dnekclock()-etime1)
#endif
      return
      END
C
c=======================================================================
      subroutine hmholtza(name,u,rhs,h1,h2,mask,mult,imsh,tli,maxit,isd)
      INCLUDE 'SIZE'
      INCLUDE 'CTIMER'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      include 'FDMH1'
      include 'TSTEP'
      include 'ADAPT'

      CHARACTER      NAME*4
      REAL           U    (LX1U*LY1U*LZ1U)
      REAL           RHS  (LX1U*LY1U*LZ1U)
      REAL           H1   (LX1U*LY1U*LZ1U)
      REAL           H2   (LX1U*LY1U*LZ1U)
      REAL           MASK (LX1U*LY1U*LZ1U)
      REAL           MULT (LX1U*LY1U*LZ1U)

      logical iffdm
      character*3 nam3
      integer  ptr,ptr1,ptr2

      tol = abs(tli)


      if (icalld.eq.0) thmhz=0.0

      iffdm = .false.
c     iffdm = .true.
      if (ifsplit) iffdm = .true.
c
      if (icalld.eq.0.and.iffdm) call set_fdm_prec_h1A

c
      icalld=icalld+1
      nhmhz=icalld
      etime1=dnekclock()
      ntot = ntota(ifield)
C     Determine which field is being computed for FDM based preconditioner bc's
c
      call chcopy(nam3,name,3)
c
                          kfldfdm = -1
c     if (nam3.eq.'TEM' ) kfldfdm =  0
c     if (name.eq.'TEM1') kfldfdm =  0  ! hardcode for temp only, for mpaul
c     if (name.eq.'VELX') kfldfdm =  1
c     if (name.eq.'VELY') kfldfdm =  2
c     if (name.eq.'VELZ') kfldfdm =  3
      if (name.eq.'PRES') kfldfdm =  ndim+1
c     if (.not.iffdm) kfldfdm=-1
C
      call dssum(rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,ntot)
      if (nid.eq.0.and.istep.le.10) 
     $    write(6,*) param(22),' p22 ',istep,imsh
      if (param(22).eq.0.or.istep.le.10)
     $    call chktcg1a(tol,rhs,h1,h2,mask,mult,imsh,isd)

      if (tli.lt.0) tol=tli ! caller-specified relative tolerance

      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr+ntota(1)
      if (imsh.eq.1) call cggoa
     $   (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1a,name)
      if (imsh.eq.2) call cggoa
     $   (u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1a(ptr),name)


      thmhz=thmhz+(dnekclock()-etime1)
      return
      END
C
c=======================================================================
      subroutine chktcg1a(tol,res,h1,h2,mask,mult,imesh,isd)
C-------------------------------------------------------------------
C
C     Check that the tolerances are not too small for the CG-solver.
C     Important when calling the CG-solver (Gauss-Lobatto mesh) with
C     zero Neumann b.c.
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'EIGEN'
      INCLUDE 'ADAPT'
      
      COMMON /ISTEP2/ IFIELD 
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      COMMON /CTMP0/ W1   (LX1U*LY1U*LZ1U*LELT)
     $ ,             W2   (LX1U*LY1U*LZ1U*LELT)
      REAL RES  (LX1U*LY1U*LZ1U)
      REAL H1   (LX1U*LY1U*LZ1U)
      REAL H2   (LX1U*LY1U*LZ1U)
      REAL MULT (LX1U*LY1U*LZ1U)
      REAL MASK (LX1U*LY1U*LZ1U)
      integer ptr,ptr1,ptr2 
C
      IF (EIGAA.NE.0.) THEN
         ACONDNO = EIGGA/EIGAA
      ELSE
         ACONDNO = 10.
      ENDIF

C
C     Single or double precision???
C
      DELTA = 1.E-9
      X     = 1.+DELTA
      Y     = 1.
      DIFF  = ABS(X-Y)
      IF (DIFF.EQ.0.) EPS = 1.E-6
      IF (DIFF.GT.0.) EPS = 1.E-13
C
      IF (IMESH.EQ.1) THEN
          NL  = NELV
          VOL = VOLVM1
      ELSEIF (IMESH.EQ.2) THEN
          NL  = NELT
          VOL = VOLTM1
      ENDIF
      NTOT1 = ntota(ifield)
      CALL COPY (W1,RES,NTOT1)
C
      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr + ntota(1)
      ptr1 = adptr(1,ifield)
      IF (IMESH.EQ.1) THEN
         CALL COL3 (W2,BINVM1A(ptr),W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT(ptr1),NTOT1)/VOLVM1)
      ELSE
         CALL COL3 (W2,BINVM1A(ptr),W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT(ptr1),NTOT1)/VOLTM1)
      ENDIF
      RMIN   = EPS*RINIT
      IF (TOL.LT.RMIN) THEN
         IF (NID.EQ.0.AND.IFPRINT)
     $   WRITE (6,*) 'New CG1-tolerance (RINIT*epsm) = ',RMIN,TOL
         TOL = RMIN
      ENDIF
C
      ptr1 = adptr(1,ifield)
      CALL RONE (W1,NTOT1)
      BCNEU1 = GLSC3(W1,MASK(ptr1),MULT(ptr1),NTOT1)
      BCNEU2 = GLSC3(W1,W1  ,MULT(ptr1),NTOT1)
      BCTEST = ABS(BCNEU1-BCNEU2)
C
      CALL AXHELMA(W2,W1,H1,H2,IMESH,ISD)
      CALL COL2   (W2,W2,NTOT1)
      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr + ntota(1)
      CALL COL2   (W2,BM1A(ptr),NTOT1)
      BCROB  = SQRT(GLSUM(W2,NTOT1)/VOL)
C
      IF ((BCTEST .LT. .1).AND.(BCROB.LT.(EPS*ACONDNO))) THEN
C         OTR = GLSC3 (W1,RES,MULT,NTOT1)
         TOLMIN = RINIT*EPS*10.
         IF (TOL .LT. TOLMIN) THEN
             TOL = TOLMIN
             IF (NID.EQ.0.AND.IFPRINT)
     $       WRITE(6,*) 'New CG1-tolerance (Neumann) = ',TOLMIN
         ENDIF
      ENDIF
C
      return
      end
c=======================================================================
      subroutine cggoa(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name)
C-------------------------------------------------------------------------
C
C     Solve the Helmholtz equation, H*U = RHS,
C     using preconditioned conjugate gradient iteration.
C     Preconditioner: diag(H).
C
C------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      include 'FDMH1'
      include 'GEOM'
      include 'ADAPT'
c
      COMMON  /CPRINT/ IFPRINT, IFHZPC
      LOGICAL          IFPRINT, IFHZPC
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      logical ifmcor,ifprint_hmh
C
      real x(1),f(1),h1(1),h2(1),mask(1),mult(1),binv(1)
      parameter        (lg=lx1u*ly1u*lz1u*lelt)
      COMMON /SCRCG/ d (lg) , scalar(2)
      common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)
c
      parameter (maxcg=900)
      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /iterhm/ niterhm
      character*4 name
      integer ptr,ptr1,ptr2 
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
c
      if (ifsplit.and.name.eq.'PRES'.and.param(42).eq.0) then
         n = nx1*ny1*nz1*nelv
         call copy      (x,f,n)
         call hmh_gmres (x,h1,h2,mult,iter)
         niterhm = iter
         return
      endif
c      write(6,*) ifsplit,name,param(44),' P44 C'

c **  zero out stuff for Lanczos eigenvalue estimator
      call rzero(diagt,maxcg)
      call rzero(upper,maxcg)
C
C     Initialization
C
      NEL    = NELV
      VOL    = VOLVM1
      IF (IMSH.EQ.2) NEL=NELT
      IF (IMSH.EQ.2) VOL=VOLTM1
      n = ntota(ifield) 
c
      tol=abs(tin)
      if (param(22).ne.0) tol=abs(param(22))
      if (name.eq.'PRES'.and.param(21).ne.0) tol=abs(param(21))
      if (tin.lt.0)       tol=abs(tin)
      niter = min(maxit,maxcg)

C     Speed-up for undeformed elements and constant properties.
      if (.not.ifsolv) then
         call setfast(h1,h2,imesh)
         ifsolv = .true.
      endif
C
C     Set up diag preconditioner.
C
      if (kfldfdm.lt.0) then
         call setpreca(D,h1,h2,imsh,isd)
      elseif(param(100).ne.2) then
         call set_fdm_prec_h1b(d,h1,h2,nel)
      endif
      
      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

!      if (nid==1) then
!        print *, 'D'
!        print '(25F9.4)', d
!        print *
!      end if
c
c     Check for non-trivial null-space
c
      ifmcor = .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)
      if (skmin.gt.0.and.h2max.eq.0) ifmcor = .true.
C
      if (name.eq.'PRES') then
c        call ortho (r)           ! Commented out March 15, 2011,pff
      elseif (ifmcor) then
         smean = -1./glsum(mult,n)
         rmean = smean*glsc2(r,mult,n)
         call cadd(r,rmean,n)
      endif
C
      krylov = 0
      rtz1=1.0
      niterhm = 0

      do iter=1,niter
Ci        
         if (kfldfdm.lt.0) then  ! Jacobi Preconditioner                      ! Step 7
c           call copy(z,r,n)
            call col3(z,r,d,n)
         else                                       ! Schwarz Preconditioner
            if (name.eq.'PRES'.and.param(100).eq.2) then
               call h1_overlap_2(z,r,mask)
               call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
               call add2         (z,w,n)
            else   
               call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
               if (name.eq.'PRES') then 
                 call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
                 call add2         (z,w,n)
               endif
            endif
         endif
!         if (nid==1) then
!           print *, 'Z'
!           print '(25F9.4)', z
!           print *
!         end if
c
         if (name.eq.'PRES') then
            call ortho (z)
         elseif (ifmcor) then
            rmean = smean*glsc2(z,mult,n)
            call cadd(z,rmean,n)
         endif
c
         rtz2=rtz1
         scalar(1)=vlsc3 (z,r,mult,n)
         scalar(2)=vlsc32(r,mult,binv,n)
         call gop(scalar,w,'+  ',2)
         rtz1=scalar(1)
         rbn2=sqrt(scalar(2)/vol)

!         if (nid==1) then
!           print *, 'SCALARS'
!           print *, scalar
!           print *
!         end if        
 
         if (iter.eq.1) rbn0 = rbn2
         if (param(22).lt.0) tol=abs(param(22))*rbn0
         if (tin.lt.0)       tol=abs(tin)*rbn0

         ifprint_hmh = .false.
         if (nid.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
         if (nid.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.

         if (ifprint_hmh)
     $      write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1),h2(1)


c        Always take at least one iteration   (for projection) pff 11/23/98
#ifndef TST_WSCAL
         IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
#else
         iter_max = param(150)
         if (name.eq.'PRES') iter_max = param(151)
         if (iter.gt.iter_max) then
#endif
c        IF (rbn2.LE.TOL) THEN
            NITER = ITER-1
c           IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
            if (nid.eq.0)
     $         write(6,3000) istep,name,niter,rbn2,rbn0,tol
            goto 9999
         ENDIF
c
         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0                 ! 8
         call add2s1 (p,z,beta,n)                ! 9

!         if (nid==1) then
!           print *, 'rtz1, rtz2, beta'
!           print *, rtz1, rtz2, beta
!           print *
!         end if

         debug = .true.
         call axhelma(w,p,h1,h2,imsh,isd)        ! 3
!         if (nid==1) then
!           print *, 'W'
!           print '(25F9.4)', w
!           print *
!         end if        
         call dssum(w,nx1,ny1,nz1)             ! 3
         call col2   (w,mask,n)                  ! 3
!         if (nid==1) then
!           print *, 'W again'
!           print '(25F9.4)', w
!           print *
!         end if        
c
         rho0 = rho
!         if (nid==0) then
!           print *, 'NID', nid
!           print *
!           print *, 'W'
!           print '(16F9.4)', w
!           print *
!           print *, 'P'
!           print '(16F9.4)', p
!           print *
!           print *, 'MULT'
!           print '(16F9.4)', mult
!           print *
!           print *, 'N', n
!           print *
!         end if
!         call mpi_barrier(nekcomm,ierr)
!         if (nid==1) then
!           print *, 'NID', nid
!           print *
!           print *, 'W'
!           print '(16F9.4)', w
!           print *
!           print *, 'P'
!           print '(16F9.4)', p
!           print *
!           print *, 'MULT'
!           print '(16F9.4)', mult
!           print *
!           print *, 'N', n
!           print *
!         end if
!         call mpi_barrier(nekcomm,ierr)
!         if (nid==2) then
!           print *, 'NID', nid
!           print *
!           print *, 'W'
!           print '(16F9.4)', w
!           print *
!           print *, 'P'
!           print '(16F9.4)', p
!           print *
!           print *, 'MULT'
!           print '(16F9.4)', mult
!           print *
!           print *, 'N', n
!           print *
!         end if
!         call mpi_barrier(nekcomm,ierr)
         rho  = glsc3(w,p,mult,n)                ! 4
         alpha=rtz1/rho                          ! 4
         alphm=-alpha
!         if (nid==1) then
!           print *, 'ALPHA'
!           print *, alpha
!           print *
!         end if        
         call add2s2(x,p ,alpha,n)               ! 5
         call add2s2(r,w ,alphm,n)               ! 6
c
c        Generate tridiagonal matrix for Lanczos scheme
         if (iter.eq.1) then
            krylov = krylov+1
            diagt(iter) = rho/rtz1
         elseif (iter.le.maxcg) then
            krylov = krylov+1
            diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
            upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
         endif
      enddo
      niter = iter-1
c
      if (nid.eq.0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
 3000 format(4x,i7,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
 3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
 3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
 9999 continue
      niterhm = niter
      ifsolv = .false.
c
c
c     Call eigenvalue routine for Lanczos scheme:
c          two work arrays are req'd if you want to save "diag & upper"
c
c     if (iter.ge.3) then
c        niter = iter-1
c        call calc (diagt,upper,w,z,krylov,dmax,dmin)
c        cond = dmax/dmin
c        if (nid.eq.0) write(6,6) istep,cond,dmin,dmax,' lambda'
c     endif
c   6 format(i9,1p3e12.4,4x,a7)
c
c     if (n.gt.0) write(6,*) 'quit in cggo'
c     if (n.gt.0) call exitt
c     call exitt
      return
      end
c=======================================================================
      subroutine setpreca(dpcm1,helm1,helm2,imesh,isd)
C-------------------------------------------------------------------
C
C     Generate diagonal preconditioner for the Helmholtz operator.
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'ADAPT'

      REAL            DPCM1 (LX1U*LY1U*LZ1U)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV, GEOFLAG
      COMMON /ISTEP2/ IFIELD
      REAL            HELM1(LX1U*LY1U*LZ1U), HELM2(LX1U*LY1U*LZ1U)
      REAL YSM1(LY1)
      INTEGER nx1o, ind,ind1,ind2,ind3,ind4,ind5,ptrvec(6),ctr
      integer ptr,ptr1,ptr2 

      IF(IMESH.EQ.1) NEL=NELV
      IF(IMESH.EQ.2) NEL=NELT
      
      NTOT = ntota(ifield)
c     The following lines provide a convenient debugging option
c     call rone(dpcm1,ntot)
c     return

      CALL RZERO(DPCM1,NTOT)
      DO 1000 IE=1,NEL
         nx1o = nx1;
         if (ie.eq.1) nx1o = 0
         call getord(ie,ifield)
         ! compute new derivative matrices if needed
         geoflag = .false.
         if (nx1.ne.lx1) then    
           geoflag = .true. 
           call ComputeG(ie,GM1A)
         else
           if (nx1o.eq.nx1) call ComputeDxyz 
         end if

        IF (IFAXIS) CALL SETAXDY ( IFRZER(IE) )
 
        NXY = NX1*NY1
        NXYZ = NX1*NY1*NZ1
        ptr = adptr1(ie,ifield)
        do iq=1,6
          ptrvec(iq) = (iq-1)*nxyz+1
        end do 
        
        if (geoflag) then  
          DO 320 IQ=1,NX1
          ctr = 0
          DO 320 IZ=1,NZ1
          DO 320 IY=1,NY1
          ind2 = ptrvec(1)+(IZ-1)*NXY+(IY-1)*NX1+(IQ-1)
          DO 320 IX=1,NX1
             ind1 = (IQ-1)*NX1+IX
             DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                            GM1A(ind2) * DXTM1E(ind1)**2
              ctr = ctr + 1
  320        CONTINUE
          DO 340 IQ=1,NY1
          ctr = 0
          DO 340 IZ=1,NZ1
          DO 340 IY=1,NY1
          ind1 = (IQ-1)*NY1+IY
          DO 340 IX=1,NX1
             ind2 = ptrvec(2)+(IZ-1)*NXY+(IQ-1)*NX1+(IX-1)
             DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                            GM1A(ind2) * DYTM1E(ind1)**2
              ctr = ctr + 1
  340        CONTINUE
          IF (NDIM.EQ.3) THEN
             DO 360 IQ=1,NZ1
             ctr = 0
             DO 360 IZ=1,NZ1
             ind1 = (IQ-1)*NZ1+IZ
             DO 360 IY=1,NY1
             DO 360 IX=1,NX1
                ind2 = ptrvec(3)+(IQ-1)*NXY+(IY-1)*NX1+(IX-1)
                DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                               GM1A(ind2) * DZTM1E(ind1)**2
                ctr = ctr + 1
  360        CONTINUE
C
C         Add cross terms if element is deformed.
C
          IF (IFDFRM(IE)) THEN
             NXX = NX1*NX1
             DO 600 IY=1,NY1
             ind2 = (IY-1)*NY1+IY
             DO 600 IZ=1,NZ1
             ind3 = (IZ-1)*NZ1+IZ
             ind1 = ptr + (IZ-1)*NXY+(IY-1)*NX1
             ind4 = ind1-ptr+ptrvec(4)
             ind5 = ind1-ptr+ptrvec(5)
             DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DXTM1E(1)*DYTM1E(ind2)
     $              + GM1A(ind5) * DXTM1E(1)*DZTM1E(ind3)
             ind1 = ind1+(NX1-1)
             ind4 = ind4+(NX1-1)
             ind5 = ind5+(NX1-1)
             DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DXTM1E(NXX)*DYTM1E(ind2)
     $              + GM1A(ind5) * DXTM1E(NXX)*DZTM1E(ind3)
  600        CONTINUE
             NXX = NY1*NY1
             DO 700 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 700 IZ=1,NZ1
               ind3 = (IZ-1)*NZ1+IZ
               ind1 = ptr + (IZ-1)*NXY+(IX-1)
               ind4 = ind1-ptr+ptrvec(4)
               ind5 = ind1-ptr+ptrvec(6)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DYTM1E(1)*DXTM1E(ind2)
     $              + GM1A(ind5) * DYTM1E(1)*DZTM1E(ind3)
               ind1 = ptr + (IZ-1)*NXY+(NY1-1)*NX1+(IX-1)
               ind4 = ind1-ptr+ptrvec(4)
               ind5 = ind1-ptr+ptrvec(6)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DYTM1E(NXX)*DXTM1E(ind2)
     $              + GM1A(ind5) * DYTM1E(NXX)*DZTM1E(ind3)
  700        CONTINUE
             NXX = NZ1*NZ1
             DO 800 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 800 IY=1,NY1
               ind3 = (IY-1)*NY1+IY
               ind1 = ptr + (IY-1)*NX1+(IX-1)
               ind4 = ind1-ptr+ptrvec(5)
               ind5 = ind1-ptr+ptrvec(6)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DZTM1E(1)*DXTM1E(ind2)
     $              + GM1A(ind5) * DZTM1E(1)*DYTM1E(ind3)
               ind1 = ptr + (NZ1-1)*NXY+(IY-1)*NX1+(IX-1)
               ind4 = ind1-ptr+ptrvec(5)
               ind5 = ind1-ptr+ptrvec(6)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind4) * DZTM1E(NXX)*DXTM1E(ind2)
     $              + GM1A(ind5) * DZTM1E(NXX)*DYTM1E(ind3)
  800        CONTINUE
          ENDIF
        ELSE
C
         IF (IFDFRM(IE)) THEN
             IZ=1
             NXX = NX1*NX1
             DO 602 IY=1,NY1
               ind2 = (IY-1)*NY1+IY
               ind1 = ptr + (IY-1)*NX1+(IX-1)
               ind3 = ind1-ptr+ptrvec(4)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind3) * DXTM1E(1)*DYTM1E(ind2)
               ind1 = ptr + (IY-1)*NX1+(NX1-1)
               ind3 = ind1-ptr+ptrvec(4)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind3) * DXTM1E(NXX)*DYTM1E(ind2)
  602        CONTINUE
             NXX = NY1*NY1
             DO 702 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 702 IZ=1,NZ1
               ind1 = ptr + (IZ-1)*NXY+(IX-1)
               ind3 = ind1-ptr+ptrvec(4)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind3) * DYTM1E(1)*DXTM1E(ind2)
               ind1 = ptr + (IZ-1)*NXY+(NY1-1)*NX1+(IX-1)
               ind3 = ind1-ptr+ptrvec(4)
               DPCM1(ind1) = DPCM1(ind1)
     $              + GM1A(ind3) * DYTM1E(NXX)*DXTM1E(ind2)
  702        CONTINUE
           ENDIF
        ENDIF
C
  
        else    ! Geoflag = .false.

          DO 520 IQ=1,NX1
          ctr = 0
          DO 520 IZ=1,NZ1
          DO 520 IY=1,NY1
          DO 520 IX=1,NX1
             ind1 = (IQ-1)*NX1+IX
             DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                            G1M1(IQ,IY,IZ,IE) * DXTM1E(ind1)**2
              ctr = ctr + 1
  520        CONTINUE
          DO 540 IQ=1,NY1
          ctr = 0
          DO 540 IZ=1,NZ1
          DO 540 IY=1,NY1
          ind1 = (IQ-1)*NY1+IY
          DO 540 IX=1,NX1
             DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                            G2M1(IX,IQ,IZ,IE) * DYTM1E(ind1)**2
              ctr = ctr + 1
  540        CONTINUE
          IF (NDIM.EQ.3) THEN
             DO 560 IQ=1,NZ1
             ctr = 0
             DO 560 IZ=1,NZ1
             ind1 = (IQ-1)*NZ1+IZ
             DO 560 IY=1,NY1
             DO 560 IX=1,NX1
                DPCM1(ptr+ctr) = DPCM1(ptr+ctr) + 
     $                               G3M1(IX,IY,IQ,IE) * DZTM1E(ind1)**2
                ctr = ctr + 1
  560        CONTINUE
C
C         Add cross terms if element is deformed.
C
          IF (IFDFRM(IE)) THEN
             NXX = NX1*NX1
             DO 620 IY=1,NY1
             ind2 = (IY-1)*NY1+IY
             DO 620 IZ=1,NZ1
             ind3 = (IZ-1)*NZ1+IZ
             ind1 = ptr + (IZ-1)*NXY+(IY-1)*NX1
             DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(1,IY,IZ,IE) * DXTM1E(1)*DYTM1E(ind2)
     $              + G5M1(1,IY,IZ,IE) * DXTM1E(1)*DZTM1E(ind3)
             ind1 = ind1+(NX1-1)
             DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(NX1,IY,IZ,IE) * DXTM1E(NXX)*DYTM1E(ind2)
     $              + G5M1(NX1,IY,IZ,IE) * DXTM1E(NXX)*DZTM1E(ind3)
  620        CONTINUE
             NXX = NY1*NY1
             DO 720 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 720 IZ=1,NZ1
               ind3 = (IZ-1)*NZ1+IZ
               ind1 = ptr + (IZ-1)*NXY+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(IX,1,IZ,IE) * DYTM1E(1)*DXTM1E(ind2)
     $              + G6M1(IX,1,IZ,IE) * DYTM1E(1)*DZTM1E(ind3)
               ind1 = ptr + (IZ-1)*NXY+(NY1-1)*NX1+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(IX,NY1,IZ,IE) * DYTM1E(NXX)*DXTM1E(ind2)
     $              + G6M1(IX,NY1,IZ,IE) * DYTM1E(NXX)*DZTM1E(ind3)
  720        CONTINUE
             NXX = NZ1*NZ1
             DO 820 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 820 IY=1,NY1
               ind3 = (IY-1)*NY1+IY
               ind1 = ptr + (IY-1)*NX1+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G5M1(IX,IY,1,IE) * DZTM1E(1)*DXTM1E(ind2)
     $              + G6M1(IX,IY,1,IE) * DZTM1E(1)*DYTM1E(ind3)
               ind1 = ptr + (NZ1-1)*NXY+(IY-1)*NX1+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G5M1(IX,IY,NZ1,IE) * DZTM1E(NXX)*DXTM1E(ind2)
     $              + G6M1(IX,IY,NZ1,IE) * DZTM1E(NXX)*DYTM1E(ind3)
  820        CONTINUE
          ENDIF
        ELSE
C
         IF (IFDFRM(IE)) THEN
             IZ=1
             NXX = NX1*NX1
             DO 608 IY=1,NY1
               ind2 = (IY-1)*NY1+IY
               ind1 = ptr + (IY-1)*NX1+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(1,IY,IZ,IE) * DXTM1E(1)*DYTM1E(ind2)
               ind1 = ptr + (IY-1)*NX1+(NX1-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(NX1,IY,IZ,IE) * DXTM1E(NXX)*DYTM1E(ind2)
  608        CONTINUE
             NXX = NY1*NY1
             DO 708 IX=1,NX1
             ind2 = (IX-1)*NX1+IX
             DO 708 IZ=1,NZ1
               ind1 = ptr + (IZ-1)*NXY+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(IX,1,IZ,IE) * DYTM1E(1)*DXTM1E(ind2)
               ind1 = ptr + (IZ-1)*NXY+(NY1-1)*NX1+(IX-1)
               DPCM1(ind1) = DPCM1(ind1)
     $              + G4M1(IX,NY1,IZ,IE) * DYTM1E(NXX)*DXTM1E(ind2)
  708        CONTINUE
           ENDIF
        ENDIF

      end if  ! End if geoflag

 1000   CONTINUE
C

      CALL COL2    (DPCM1,HELM1,NTOT)
      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr+ntota(1)
      CALL ADDCOL3 (DPCM1,HELM2,BM1A(ptr),NTOT)
C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
C      IF (IFAXIS.AND.(ISD.EQ.2)) THEN
C         DO 1200 IEL=1,NEL
CC
C            IF (IFRZER(IEL)) THEN
C               CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
C            ENDIF
CC
C            DO 1190 J=1,NY1
C            DO 1190 I=1,NX1
C               IF (YM1(I,J,1,IEL).NE.0.) THEN
C                  TERM1 = BM1(I,J,1,IEL)/YM1(I,J,1,IEL)**2
C                  IF (IFRZER(IEL)) THEN
C                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)
C     $                       *JACM1(I,1,1,IEL)/YSM1(I)
C                  ELSE
C                     TERM2 = 0.
C                  ENDIF
C                  DPCM1(I,J,1,IEL) = DPCM1(I,J,1,IEL)
C     $                             + HELM1(I,J,1,IEL)*(TERM1+TERM2)
C               ENDIF
C 1190       CONTINUE
C 1200    CONTINUE
C      ENDIF
C
      CALL DSSUM (DPCM1,NX1,NY1,NZ1)
      CALL INVCOL1 (DPCM1,NTOT)
C
      return
      END
C
c------------------------------------------------------------------
      subroutine sethlma(h1,h2,intloc)
 
c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTLOC =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'ADAPT'

      real h1(1),h2(1)
      integer ptr,ptr1,ptr2 

      nel   = nelfld(ifield)
      ntot1 = ntota(ifield) 
      ptr = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr+ntota(1)

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiffa(ptr),ntot1)
         if (intloc.eq.0) then
            call rzero (h2,ntot1)
         else
            if (ifield.eq.1.or.param(107).eq.0) then 
               call cmult2 (h2,vtransa(ptr),dtbd,ntot1)
            else   ! unsteady reaction-diffusion type equation

               do i=1,ntot1
                 h2(i) = dtbd*vtransa(ptr+i) + param(107)
               enddo

            endif

         endif

c        if (ifield.eq.1 .and. ifanls) then   ! this should be replaced
c           const = 2.                        ! with a correct stress
c           call cmult (h1,const,ntot1)       ! formulation
c        endif

      ELSE
         CALL COPY  (H1,VDIFFA(ptr),NTOT1)
         CALL RZERO (H2,NTOT1)
         if (param(107).ne.0) then
            write(6,*) 'SPECIAL SETHLM!!',param(107)
c           call cfill (h2,param(107),ntot1)
            call copy  (h2,vtransa(ptr),ntot1)
         endif
      ENDIF

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE BCNEUSCA(S,ITYPE)
C
C     Apply Neumann boundary conditions to surface of scalar, S.
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
C     If ITYPE = 1, then S is returned as the rhs contribution to the 
C                   volumetric flux.
C
C     If ITYPE =-1, then S is returned as the lhs contribution to the 
C                   diagonal of A.
C
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'CTIMER'
      INCLUDE 'NEKUSE'
C
      DIMENSION S(LX1U*LY1U*LZ1U*LELT)
      common  /nekcb/ cb
      CHARACTER CB*3
      integer bptr1,ctr,ix,iy,iz
      integer ptr,ptr1,ptr2,ptr3 
C
#ifndef NOTIMER
      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dnekclock()
#endif
C
      NFACES=2*NDIM
      NXYZ  =NX1*NY1*NZ1
      NEL   =NELFLD(IFIELD)
      ntot = ntota(ifield)
      CALL RZERO(S,NTOT)
C
      IF (ITYPE.EQ.-1) THEN
C
C        Compute diagonal contributions to accomodate Robin boundary conditions
C
         DO 1000 IE=1,NEL
          call getord(ie,ifield)
          NXYZ  =NX1*NY1*NZ1
          NXY  = NX1*NY1
          NXZ  = NX1*NZ1
          NYZ  = NY1*NZ1
          ptr1 = adptr(ie,ifield)-1
          ptr2 = adptr1(ie,ifield)-1
          ptr3 = ptr1
          if (ifield.ne.1) ptr3 = ptr3+ntota(1) ! -1
         DO 1000 IFACE=1,NFACES
            ieg=lglel(ie)
            CB =CBC(IFACE,IE,IFIELD)
            IF (CB.EQ.'C  ' .OR. CB.EQ.'c  ' .OR.
     $          CB.EQ.'R  ' .OR. CB.EQ.'r  ') THEN
C
               IF (CB.EQ.'C  ') HC   = BC(2,IFACE,IE,IFIELD)
               IF (CB.EQ.'R  ') THEN
                                TINF = BC(1,IFACE,IE,IFIELD)
                                HRAD = BC(2,IFACE,IE,IFIELD)
               ENDIF
               IA=0
C
C IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
C
               CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
               ! NEED TO CHANGE THIS SECTION FOR ADAPTIVE CASE
               DO 100 IZ=KZ1,KZ2
               DO 100 IY=KY1,KY2
               DO 100 IX=KX1,KX2
                  IA = IA + 1
                  TS = T(ptr1+IA,1,1,1,1)
                  IF (CB.EQ.'c  ' .OR. CB.EQ.'r  ') THEN
                     CALL NEKASGN (ptr2+IA,1,1,1) 
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'r  ' .OR. CB.EQ.'R  ') 
     $               HC = HRAD * (TINF**2 + TS**2) * (TINF + TS)
C                  S(ptr2+IA) = S(ptr2+IA) +
C     $               HC*AREA(IA,1,IFACE,IE)/BM1A(ptr3+IA)
                  S(ptr2+IA) = S(ptr2+IA) +
     $               HC*AREAA(bptr1+IA+(IFACE-1)*nxz-1)/BM1A(ptr3+IA)   ! No area for now
  100          CONTINUE
            ENDIF
 1000    CONTINUE
      ENDIF
      IF (ITYPE.EQ.1) THEN
C
C        Add passive scalar fluxes to rhs
C
         DO 2000 IE=1,NEL
          call getord(ie,ifield)
          ptr1 = adptr(ie,ifield)-1
          ptr2 = adptr1(ie,ifield)-1
          bptr1 = bptr(ie,ifield)
          nxz = nx1*nz1
          nxy = nx1*ny1
          nyz = ny1*nz1
          nxyz = nx1*ny1*nz1
         DO 2000 IFACE=1,NFACES
            ieg=lglel(ie)
            CB =CBC(IFACE,IE,IFIELD)
            IF (CB.EQ.'F  ' .OR. CB.EQ.'f  ' .OR.
     $          CB.EQ.'C  ' .OR. CB.EQ.'c  ' .OR. 
     $          CB.EQ.'R  ' .OR. CB.EQ.'r  ' ) THEN
C
                IF (CB.EQ.'F  ') FLUX=BC(1,IFACE,IE,IFIELD)
                IF (CB.EQ.'C  ') FLUX=BC(1,IFACE,IE,IFIELD)
     $                               *BC(2,IFACE,IE,IFIELD)
                IF (CB.EQ.'R  ') THEN
                                 TINF=BC(1,IFACE,IE,IFIELD)
                                 HRAD=BC(2,IFACE,IE,IFIELD)
                ENDIF
C
C              Add local weighted flux values to rhs, S.
C
C IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
               IA=0
               CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
               DO 200 IZ=KZ1,KZ2
               DO 200 IY=KY1,KY2
               DO 200 IX=KX1,KX2
                  IA = IA + 1
                  ctr = (IZ-1)*NXY+(IY-1)*NX1+IX 
                  !TS = Tad(ptr1+ctr)
                  IF (CB.EQ.'f  ') THEN
                     CALL NEKASGN (ctr,1,1,ie)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'c  ') THEN
                     CALL NEKASGN (ctr,1,1,ie)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                     FLUX = TINF*HC
                  ENDIF
                  IF (CB.EQ.'r  ') THEN
                     CALL NEKASGN (ctr,1,1,ie)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'R  ' .OR. CB.EQ.'r  ') 
     $               FLUX = HRAD*(TINF**2 + TS**2)*(TINF + TS) * TINF
C
C                 Add computed fluxes to boundary surfaces:
                  S(ptr2+ctr) = S(ptr2+ctr)
     $                 + FLUX*AREAA(bptr1+ia+(IFACE-1)*nxz-1)
  200          CONTINUE
            ENDIF
 2000    CONTINUE
      ENDIF

C
#ifndef NOTIMER
      tusbc=tusbc+(dnekclock()-etime1)
#endif
C
      RETURN
      END
c------------------------------------------------------------------
      SUBROUTINE BCDIRSCA(S)
C
C     Apply Dirichlet boundary conditions to surface of scalar, S.
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
      INCLUDE 'CTIMER'
      INCLUDE 'ADAPT'
C
      DIMENSION S(LX1U*LY1U*LZ1U*LELT)
      COMMON /SCRSFA/ TMP(LX1U*LY1U*LZ1U*LELT,1,1,1)
     $             , TMA(LX1U*LY1U*LZ1U*LELT)
     $             , SMU(LX1U*LY1U*LZ1U*LELT)
      common  /nekcb/ cb
      CHARACTER CB*3
      integer bcptr
      integer ptr,ptr1,ptr2 
 
#ifndef NOTIMER
      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dnekclock()
#endif
C
      IFLD   = 1
      NFACES = 2*NDIM
      NEL    = NELFLD(IFIELD)
      NFLDT  = NFIELD - 1
      ntot = ntota(ifield)
C
      CALL RZERO(TMP,NTOT)
C
C     Temperature boundary condition
C
      DO 2100 ISWEEP=1,2
C
         IF (IFMODEL .AND. IFKEPS .AND. IFIELD.GE.NFLDT)
     $       CALL TURBWBC (TMP,TMA,SMU)  
C
         DO 2010 IE=1,NEL
         call getord(ie,ifield)
         bcptr = adptr1(ie,ifield)
         iev = 1
         DO 2010 IFACE=1,NFACES
            CB=CBC(IFACE,IE,IFIELD)
            BC1=BC(1,IFACE,IE,IFIELD)
            BC2=BC(2,IFACE,IE,IFIELD)
            BC3=BC(3,IFACE,IE,IFIELD)
            BC4=BC(4,IFACE,IE,IFIELD)
            BCK=BC(4,IFACE,IE,IFLD)
            BCE=BC(5,IFACE,IE,IFLD)
            IF (CB.EQ.'T  ') CALL FACEV (TMP(bcptr,1,1,1),IEV,IFACE,BC1,
     $                                   NX1,NY1,NZ1)
            IF (CB.EQ.'MCI') CALL FACEV (TMP(bcptr,1,1,1),IEV,IFACE,BC4,
     $                                   NX1,NY1,NZ1)
            IF (CB.EQ.'MLI') CALL FACEV (TMP(bcptr,1,1,1),IEV,IFACE,BC4,
     $                                   NX1,NY1,NZ1)
            IF (CB.EQ.'KD ') CALL FACEV (TMP(bcptr,1,1,1),IEV,IFACE,BCK,
     $                                   NX1,NY1,NZ1)
            IF (CB.EQ.'ED ') CALL FACEV (TMP(bcptr,1,1,1),IEV,IFACE,BCE,
     $                                   NX1,NY1,NZ1)
            IF (CB.EQ.'t  ' .OR. CB.EQ.'kd ' .OR. CB.EQ.'ed ') then
               CALL FACEISA(CB,TMP(bcptr,1,1,1),IE,IFACE,NX1,NY1,NZ1) 
            end if
 2010    CONTINUE
C
C        Take care of Neumann-Dirichlet shared edges...
C
         IF (ISWEEP.EQ.1) CALL DSOP(TMP,'MXA',NX1,NY1,NZ1)
         IF (ISWEEP.EQ.2) CALL DSOP(TMP,'MNA',NX1,NY1,NZ1)
 2100 CONTINUE
C
C     Copy temporary array to temperature array.
C
      bcptr = adptr(1,ifield)
      CALL COL2(S,TMASKA(bcptr),NTOT)
      CALL ADD2(S,TMP,NTOT)

#ifndef NOTIMER
      tusbc=tusbc+(dnekclock()-etime1)
#endif

      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine makeuqa

c     Fill up user defined forcing function and collocate will the
c     mass matrix on the Gauss-Lobatto mesh.

      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'ADAPT'

      integer ntot
      integer ptr,ptr1,ptr2 

      ntot = ntota(ifield)
      ptr = adptr(1,ifield)
      ptr2 = ptr
      if (ifield.ne.1) ptr2 = ptr2+ntota(1)

      time = time-dt        ! Set time to t^n-1 for user function

      call rzero   ( bqa(ptr) ,    ntot)
      call setqvola( bqa(ptr)          )
      call col2    ( bqa(ptr) ,bm1a(ptr2),ntot)

      time = time+dt        ! Restore time

      return
      end
c-----------------------------------------------------------------------
      subroutine setqvola(bql)

c     Set user specified volumetric forcing function (e.g. heat source).

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'ADAPT'

      real bql(lx1u*ly1u*lz1u*lelt)
      integer lxyz1,ntot1,nxyz,iel
      integer ptr,ptr1,ptr2 

      nel   = nelfld(ifield)
      nxyz1 = nx1*ny1*nz1
      lxyz1 = lx1*ly1*lz1
      ntot1 = nxyz1*nel

      do iel=1,nel
         igrp = igroup(iel)
         if (matype(igrp,ifield).eq.1) then ! constant source within a group
            ptr = lxyz1*(iel-1)+1
            if (ifadapt) then
              call getord(iel)
              nxyz1 = nx1*ny1*nz1
              ptr = adptr(iel,ifield)
            end if
            cqvol = cpgrp(igrp,ifield,3)
            call cfill (bql(ptr),cqvol,nxyz1)
         else  !  pff 2/6/96 ............ default is to look at userq
            call getord(iel,ifield)
            nxyz1 = nx1*ny1*nz1
            ptr = adptr1(iel,ifield)
            call MapV(vxa,vya,vza,vx,vy,vz,iel,ifield)
            call nekuqa(bql(ptr),iel)
         endif
      enddo
  101 FORMAT(' Wrong material type (',I3,') for group',I3,', field',I2
     $    ,/,' Aborting in SETQVOL.')
C   
      return
      end
c-----------------------------------------------------------------------
      subroutine nekuqa(bql,iel)
C------------------------------------------------------------------
C
C     Generate user-specified volumetric source term (temp./p.s.)
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'PARALLEL'
      include 'TSTEP'
      include 'NEKUSE'
      include 'INPUT'
c
      real bql(lx1u*ly1u*lz1u)
      integer ctr
c
#ifdef MOAB
c pulling in temperature right now, since we dont have anything else
      call userq2(bql)
      return
#endif
      ielg = lglel(iel)
      call getord(iel,ifield)
      nxy = nx1*ny1
      do 10 k=1,nz1
      do 10 j=1,ny1
      do 10 i=1,nx1
         ctr = (k-1)*nxy+(j-1)*nx1+i
         call nekasgn (ctr,1,1,iel)
         qvol = 0.0
         call userq   (ctr,1,1,ielg)
         bql(ctr) = qvol
 10   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine makeqa

      include 'SIZE'
      include 'TOTAL'
      
      integer ntot
      logical  if_conv_std
      integer ptr,ptr1,ptr2 
      
      ntot = ntota(ifield)

      call makeuqa
 
         print *, 'BQ', 1, nid
         print '(25F8.3)', BQA
         print *
      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      if (ifadvc(ifield).and..not.ifchar.and.if_conv_std) call convaba
         print *, 'BQ', 2, nid
         print '(25F8.3)', BQA
         print *
      call makeabqa
         print *, 'BQ', 3, nid
         print '(25F8.3)', BQA
         print *
      call makebdqa
         print *, 'BQ', 4, nid
         print '(25F8.3)', BQA
         print *

      end
c-----------------------------------------------------------------------
      subroutine convaba
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'ADAPT'
C
      COMMON /SCRUZA/ TA(LX1U*LY1U*LZ1U*LELT)
      integer nel,ntot1
      integer ptr,ptr1,ptr2 
C

      NEL = NELFLD(IFIELD)
      NTOT1 = ntota(ifield)
      ptr = adptr(1,ifield)
      ptr1 = ptr+ntota(1)
      CALL CONVOPA (TA,Tad(ptr))
      !print *, 'TA CONV'
      !print '(25F8.3)', TA
      !print *
      ptr = adptr(1,ifield)
      ptr1 = ptr+ntota(1)
      CALL COL2    (TA,VTRANSA(ptr1),NTOT1)
      !print *, 'TA CONV 2'
      !print '(25F8.3)', TA
      !print *
      ptr = adptr(1,ifield)
      ptr1 = ptr+ntota(1)
      CALL SUBCOL3 (BQA(ptr),BM1A(ptr1),TA,NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine convopa(conv,fi)
C
C     Compute the convective term CONV for a passive scalar field FI
C     using the skew-symmetric formulation.
C     The field variable FI is defined on mesh M1 (GLL) and
C     the velocity field is assumed given.
C
C     IMPORTANT NOTE: Use the scratch-arrays carefully!!!!!
C
C     The common-block SCRNS is used in CONV1 and CONV2.
C     The common-blocks CTMP0 and CTMP1 are also used as scratch-arrays
C     since there is no direct stiffness summation or Helmholtz-solves. 
C
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
C
C     Use the common blocks CTMP0 and CTMP1 as work space.
C
C     Arrays in parameter list
C
      REAL    CONV (LX1U*LY1U*LZ1U*LELT) 
      REAL    FI   (LX1U*LY1U*LZ1U*LELT)
      integer ptr,ptr1,ptr2 

#ifndef NOTIMER
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()
#endif
C
      NTOT1 = ntota(1)
      NTOTZ = ntota(ifield)
C
      CALL RZERO  (CONV,NTOTZ)
C
      !if (param(86).ne.0.0) then  ! skew-symmetric form
      !   call convopo(conv,fi)
      !   goto 100
      !endif

c     write(6,*) istep,param(99),' CONVOP',ifpert
c     ip99 = param(99)
c     if (istep.gt.5) call exitti(' CONVOP dbg: $',ip99)
      if (param(99).eq.2.or.param(99).eq.3) then 
         call conv1da(conv,fi)  !    use dealiased form
      elseif (param(99).eq.4) then
           ptr = adptr(1,ifield)
           ptr1 = ptr+ntota(1)
           !print *, 'VXD'
           !print '(425F8.3)', vxd
           !print *
           !print *, 'VYD'
           !print '(425F8.3)', vyd
           !print *
!           call convect_new (conv,fi,.false.,vxd,vyd,vzd,.true.)
           call convect_new (conv,fi,.false.,vx,vy,vz,.false.)
           !print *, 'CONV'
           !print '(25F8.3)', conv
           !print *
           !print *, 'BM1', ptr
           !print '(25F8.3)', bm1a
           !print *
           call invcol2     (conv,bm1a(ptr1),ntotz)  ! local mass inverse
           !print *, 'CONV'
           !print '(25F8.3)', conv
           !print *
      !elseif (param(99).eq.5) then
      !   call convect_cons(conv,fi,.false.,vx,vy,vz,.false.)
      !   call invcol2     (conv,bm1,ntot1)  ! local mass inverse
      else
         call conv1a(conv,fi)  !    use the convective form
      endif

 100  continue

#ifndef NOTIMER
      tadvc=tadvc+(dnekclock()-etime1)
#endif

      return
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine conv1da(dfi,fi)
C--------------------------------------------------------------------
C
C     Compute D*FI (part of the convection operator)
C     De-aliased version 3/11/97
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      REAL           DFI (LX1U*LY1U*LZ1U*NELT) 
      REAL           FI  (LX1U*LY1U*LZ1U*NELT) 
      REAL           wk1 (LX1U*LY1U*LZ1U) 
c
      COMMON /CTMP0A/ TA1  (LX1U*LY1U*LZ1U)
     $             ,  DFID (LXDU*LYDU*LZDU) 
     $             ,  TA1D (LXDU*LYDU*LZDU) 
C
      integer icalld
      save icalld
      data icalld /0/
      logical geoflag
      integer ie, nel,nxo,ptr,ptr1,ptr2
c
      nel = nelfld(ifield)
      nxo = 0
      do ie=1,nel
        call getord(ie,ifield)
        nxy1  = nx1*ny1
        nyz1  = ny1*nz1
        nxyz1 = nx1*ny1*nz1
        ntotd = nxd*nyd*nzd
        geoflag = .false.
        if (nx1.ne.lx1) then
          geoflag = .true.
          call ComputeG(ie,GM1A) ! For now. Really don't need G, just rxm1e, rym1e, etc.
        end if
        
        ptr1 = adptr1(ie,ifield)
        ptr2 = adptr(ie,1)

        if (nx1.ne.nxo) call setmapa(nx1,nxd)

c       interpolate ta1 and vx onto larger mesh
c
        if (geoflag) then
          CALL DUDXYZA(TA1,FI(ptr1),RXM1E,SXM1E,TXM1E,DXM1E,DYTM1E,
     $                 DZTM1E,JACM1E,JACMIE)
        else
          CALL DUDXYZA(TA1,FI(ptr1),RXM1(1,1,1,IE),SXM1(1,1,1,IE),
     $                 TXM1(1,1,1,IE),DXM1,DYTM1,DZTM1,JACM1(1,1,1,IE),
     $                 JACMI(1,IE))
        end if
        call specmp(ta1d,nxd,ta1,nx1,ixadpt,ixtadpt,wk1)
        call specmp(vxd,nxd,vx(1,1,1,ie),nx1,ixadpt,ixtadpt,wk1)
        CALL COL3   (DFID,TA1D,VXD,NTOTD)

c       interpolate ta1 and vy onto larger mesh
c
        if (geoflag) then
          CALL DUDXYZA(TA1,FI(ptr1),RYM1E,SYM1E,TYM1E,DXM1E,DYTM1E,
     $                 DZTM1E,JACM1E,JACMIE)
        else
          CALL DUDXYZA(TA1,FI(ptr1),RYM1(1,1,1,IE),SYM1(1,1,1,IE),
     $                 TYM1(1,1,1,IE),DXM1,DYTM1,DZTM1,JACM1(1,1,1,IE),
     $                 JACMI(1,IE))
        end if
        call specmp(ta1d,nyd,ta1,ny1,ixadpt,ixtadpt,wk1)
        call specmp(vyd,nyd,vy(1,1,1,ie),ny1,ixadpt,ixtadpt,wk1)
        CALL COL3   (DFID,TA1D,VYD,NTOTD)
c
        IF (if3d) THEN
c
c          interpolate ta1 and vz onto larger mesh
c
           if (geoflag) then
             CALL DUDXYZA(TA1,FI(ptr1),RZM1E,SZM1E,TZM1E,DXM1E,DYTM1E,
     $                    DZTM1E,JACM1E,JACMIE)
           else
             CALL DUDXYZA(TA1,FI(ptr1),RZM1(1,1,1,IE),SZM1(1,1,1,IE),
     $                    TZM1(1,1,1,IE),DXM1,DYTM1,DZTM1,
     $                    JACM1(1,1,1,IE),JACMI(1,IE))
           end if
           call specmp(ta1d,nzd,ta1,nz1,ixadpt,ixtadpt,wk1)
           call specmp(vzd,nzd,vz(1,1,1,ie),nz1,ixadpt,ixtadpt,wk1)
           CALL COL3   (DFID,TA1D,VZD,NTOTD)
        ENDIF
c
c       Now, *project* DFID onto mesh 1 using L2 projection
c
        if (nx1.ne.nxo) call setproj(nx1,nxd)
        call specmp(dfi(ptr1),nx1,dfid,nxd,pmd1,pmd1t,wk1)
        nxo = nx1
      end do 
      
      return
      END
C------------------------------------------------------------------------
      subroutine conv1a(du,u)
c
      include 'SIZE'
      include 'DXYZ'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'TSTEP'
      include 'ADAPT'
c
      real  du  (lx1u*ly1u*lz1u*lelt)
      real  u   (lx1u*ly1u*lz1u*lelt)
      real  vxe(lx1u*ly1u*lz1u),vye(lx1u*ly1u*lz1u),vze(lx1u*ly1u*lz1u)
c
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
C
C     Store the inverse jacobian to speed this operation up
C
      common /ctmp0a/ dudre(lx1u*ly1u*lz1u)
     $             ,  dudse(lx1u*ly1u*lz1u)
     $             ,  dudte(lx1u*ly1u*lz1u)
      logical geoflag
      integer nel,ntot,nxy1,nyz1,nxyz1
      integer ptr,ptr1,ptr2,ptr3 

      nel = nelfld(ifield)
      ntot  = ntota(ifield)
C
C     Compute vel.grad(u)
C
      do ie=1,nel
        call getord(ie,ifield)
        nxy1  = nx1*ny1
        nyz1  = ny1*nz1
        nxyz1 = nx1*ny1*nz1
        geoflag = .false.
        if (nx1.ne.lx1) then
          geoflag = .true.
          call ComputeG(ie,GM1A) ! For now. Really don't need G, just rxm1e, rym1e, etc.
          call MapV(vxe,vye,vze,vx,vy,vz,ie,ifield)  
        end if
        
        ptr1 = adptr1(ie,ifield)
        ptr2 = adptr(ie,1)
C
        !if (geoflag) then
        !  CALL DUDXYZA (dudre,U(ptr1),RXM1E,SXM1E,TXM1E,DXM1E,DYTM1E,
     $  !                DZTM1E,JACM1E,JACMIE)
        !  CALL COL3    (du(ptr1),dudre,VXE,NXYZ1)
        !  CALL DUDXYZA (dudse,U(ptr1),RYM1E,SYM1E,TYM1E,DXM1E,DYTM1E,
     $  !                DZTM1E,JACM1E,JACMIE)
        !  CALL ADDCOL3 (du(ptr1),dudse,vye,nxyz1)
        !  if (if3d) then
        !    CALL DUDXYZA (dudte,U(ptr1),RZM1E,SZM1E,TZM1E,DXM1E,DYTM1E,
     $  !                  DZTM1E,JACM1E,JACMIE)
        !    CALL ADDCOL3 (du(ptr1),dudte,VZE,NXYZ1)
        !  end if
        !else
        !  CALL DUDXYZA (dudre,U(ptr1),RXM1(1,1,1,IE),SXM1(1,1,1,IE),
     $  !                TXM1(1,1,1,IE),DXM1,DYTM1,DZTM1,JACM1(1,1,1,IE),
     $  !                JACMI(1,IE))
        !  CALL COL3    (DU(ptr1),dudre,VX(1,1,1,IE),NXYZ1)
        !  CALL DUDXYZA (dudse,U(ptr1),RYM1(1,1,1,IE),SYM1(1,1,1,IE),
     $  !                TYM1(1,1,1,IE),DXM1,DYTM1,DZTM1,JACM1(1,1,1,IE),
     $  !                JACMI(1,IE))
        !  CALL ADDCOL3 (du(ptr1),dudse,vy(1,1,1,ie),nxyz1)
        !  if (if3d) then
        !    CALL DUDXYZA (dudte,U(ptr1),RZM1(1,1,1,IE),SZM1(1,1,1,IE),
     $  !                TZM1(1,1,1,IE),DXM1,DYTM1,DZTM1,JACM1(1,1,1,IE),
     $  !                JACMI(1,IE))
        !    CALL ADDCOL3 (du(ptr1),dudte,vz(1,1,1,ie),nxyz1)
        !  end if
        !end if
c

        if (if3d) then
c
           if (geoflag) then
             call mxm   (dxm1e,nx1,u(ptr1),nx1,dudre,nyz1)
             do iz=1,nz1
               ptr3 = (iz-1)*nx1*ny1+1
               call mxm (u(ptr1+ptr3-1),nx1,dytm1e,ny1,dudse(ptr3),ny1)
             enddo
             call mxm   (u(ptr1),nxy1,dztm1e,nz1,dudte,nz1)
             do i=1,nxyz1
                du(ptr1+i-1) = jacmie(i)*(
     $                       vxe(i)*(
     $                            rxm1e(i)*dudre(i)
     $                          + sxm1e(i)*dudse(i)
     $                          + txm1e(i)*dudte(i) )
     $                     + vye(i)*(
     $                            rym1e(i)*dudre(i)
     $                          + sym1e(i)*dudse(i)
     $                          + tym1e(i)*dudte(i) )
     $                     + vze(i)*(
     $                            rzm1e(i)*dudre(i)
     $                          + szm1e(i)*dudse(i)
     $                          + tzm1e(i)*dudte(i) ) )
             enddo
           else
             call mxm   (dxm1,nx1,u(ptr1),nx1,dudre,nyz1)
             do iz=1,nz1
               ptr3 = (iz-1)*nx1*ny1+1
               call mxm (u(ptr1+ptr3-1),nx1,dytm1,ny1,dudse(ptr3),ny1)
             enddo
             call mxm   (u(ptr1),nxy1,dztm1,nz1,dudte,nz1)
             do i=1,nxyz1
                du(ptr1+i-1) = jacmi(i,ie)*(
     $                       vx(i,1,1,ie)*(
     $                            rxm1(i,1,1,ie)*dudre(i)
     $                          + sxm1(i,1,1,ie)*dudse(i)
     $                          + txm1(i,1,1,ie)*dudte(i) )
     $                     + vy(i,1,1,ie)*(
     $                            rym1(i,1,1,ie)*dudre(i)
     $                          + sym1(i,1,1,ie)*dudse(i)
     $                          + tym1(i,1,1,ie)*dudte(i) )
     $                     + vz(i,1,1,ie)*(
     $                            rzm1(i,1,1,ie)*dudre(i)
     $                          + szm1(i,1,1,ie)*dudse(i)
     $                          + tzm1(i,1,1,ie)*dudte(i) ) )
             enddo
           end if
c
        else
c
c           2D
            if (geoflag) then
              call mxm (dxm1e,nx1,u(ptr1),nx1,dudre,nyz1)
              !print *, 'ie', ie, nx1, nyz1
              !print *, 'DUDR'
              !print '(25F8.3)', dudre
              !print *
              call mxm (u(ptr1),nx1,dytm1e,ny1,dudse,ny1)
              !print *, 'DUDS'
              !print '(25F8.3)', dudse
              !print *
              !print *, 'geoflag', geoflag
              !print *
              do i=1,nxyz1
               !print '(2I8,9F9.4)',i,ie, jacmie(i),rxm1e(i),sxm1e(i),
     $         !       rym1e(i),sym1e(i),dudre(i),dudse(i),vxe(i),vye(i)
                 du(ptr1+i-1) = jacmie(i)*(
     $                        vxe(i)*(
     $                             rxm1e(i)*dudre(i)
     $                           + sxm1e(i)*dudse(i) )
     $                      + vye(i)*(
     $                             rym1e(i)*dudre(i)
     $                           + sym1e(i)*dudse(i) ) )
              enddo
            else
              call mxm (dxm1,nx1,u(ptr1),nx1,dudre,nyz1)
              !print *, 'ie', ie, nx1, nyz1
              !print *, 'dxm1'
              !print '(25F8.3)', dxm1
              !print *
              !print *, 'U', ptr1
              !print '(25F8.3)', u
              !print *
              !print *, 'IE', ie
              !print *, 'DUDR'
              !print '(25F8.3)', dudre
              !print *
              !print *, 'RMx1'
              !print '(25F8.3)', rxm1
              !print *
              !print *, 'SMx1'
              !print '(25F8.3)', sxm1
              !print *
              call mxm (u(ptr1),nx1,dytm1,ny1,dudse,ny1)
              !print *
              !print *, 'RMy1'
              !print '(25F8.3)', rym1
              !print *
              !print *, 'SMy1'
              !print '(25F8.3)', sym1
              !print *
              !print *, 'jacmi'
              !print '(25F8.3)', jacmi
              !print *
               
              !print *, 'geoflag', geoflag
              !print *
              do i=1,nxyz1
              ! print '(2I8,5F9.4)',i,ie,rxm1(i,1,1,ie),dudre(i),
     $        !      sxm1(i,1,1,ie)*dudse(i),rym1(i,1,1,ie)*dudre(i)
     $        !      +sym1(i,1,1,ie)*dudse(i),vx(i,1,1,ie),vy(i,1,1,ie),
     $        !      jacmi(i,ie)
                 du(ptr1+i-1) = jacmi(i,ie)*(
     $                        vx(i,1,1,ie)*(
     $                             rxm1(i,1,1,ie)*dudre(i)
     $                           + sxm1(i,1,1,ie)*dudse(i) )
     $                      + vy(i,1,1,ie)*(
     $                             rym1(i,1,1,ie)*dudre(i)
     $                           + sym1(i,1,1,ie)*dudse(i) ) )
              enddo
            end if
            !print *, 'DU', ptr
            !print '(25F8.3)', du
            !print *
        endif

       enddo
c
       return
       end
c-----------------------------------------------------------------------
      subroutine dudxyza(du,u,rm1,sm1,tm1,dx1,dyt1,dzt1,jm1,jmi)
C--------------------------------------------------------------
C
C     DU   - dU/dx or dU/dy or dU/dz
C     U    - a field variable defined on mesh 1
C     RM1  - dr/dx or dr/dy or dr/dz  
C     SM1  - ds/dx or ds/dy or ds/dz
C     TM1  - dt/dx or dt/dy or dt/dz
C     JM1  - the Jacobian   
C     IMESH - topology: velocity (1) or temperature (2) mesh
C
C--------------------------------------------------------------
      include 'SIZE'
      include 'DXYZ'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
C
      REAL  DU   (LX1U*LY1U*LZ1U)
      REAL  U    (LX1U*LY1U*LZ1U)
      REAL  RM1  (LX1U*LY1U*LZ1U)
      REAL  SM1  (LX1U*LY1U*LZ1U)
      REAL  TM1  (LX1U*LY1U*LZ1U)
      REAL  DX1  (LX1U*LX1U)
      REAL  DYT1 (LY1U*LY1U)
      REAL  DZT1 (LZ1U*LZ1U)
      REAL  JM1  (LX1U*LY1U*LZ1U)
      REAL  JMI  (LX1U*LY1U*LZ1U)
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV, GFLAG
C
      REAL  DRST(LX1U*LY1U*LZ1U)
      integer nxy1, nyz1, nxyz1
      integer ptr,ptr1,ptr2 
C

C
      NXY1  = NX1*NY1
      NYZ1  = NY1*NZ1
      NXYZ1 = NX1*NY1*NZ1

      IF (IFAXIS) CALL SETAXDY (IFRZER(IEL) )
C

      IF (NDIM.EQ.2) THEN
            CALL MXM     (DX1,NX1,U,NX1,DU,NYZ1)
            !print *, 'DUDR'
            !print '(25F8.3)',DU
            !print *
            CALL COL2    (DU,RM1,NXYZ1)
            !print *, 'RM1'
            !print '(25F8.3)',RM1
            !print *
            CALL MXM     (U,NX1,DYT1,NY1,DRST,NY1)
            !print *, 'DRST'
            !print '(25F8.3)',DRST
            !print *
            !print *, 'SM1'
            !print '(25F8.3)',SM1
            !print *
            !print *, 'JMI'
            !print '(25F8.3)',JMI
            !print *
            CALL ADDCOL3 (DU,DRST,SM1,NXYZ1)
      ELSE
            CALL MXM   (DX1,NX1,U,NX1,DU,NYZ1)
            CALL COL2  (DU,RM1,NXYZ1)
            DO 20 IZ=1,NZ1
               CALL MXM  (U,NX1,DYT1,NY1,DRST,NY1)
 20         CONTINUE
            CALL ADDCOL3 (DU,DRST,SM1,NXYZ1)
            CALL MXM     (U,NXY1,DZT1,NZ1,DRST,NZ1)
            CALL ADDCOL3 (DU,DRST,TM1,NXYZ1)
      ENDIF
C
 1000 CONTINUE
C
      CALL COL2 (DU,JMI,NXYZ1)
C
      return
      END
C
c-----------------------------------------------------------------------
      subroutine MapV(vxe,vye,vze,u,v,w,ie,ifld)

      include "SIZE"
      include "ADAPT"
      include 'INPUT'
      include 'SOLN'
    
      real vxe(lx1u*ly1u*lz1u),vye(lx1u*ly1u*lz1u),vze(lx1u*ly1u*lz1u),
     $     wk1(lx1u*ly1u*lz1u),wk2(lx1u*ly1u*lz1u)
      real u(lx1,ly1,lz1,1),v(lx1,ly1,lz1,1),w(lx1,ly1,lz1,1)
      integer ie,ntot,nxyz,ifld
      integer ptr,ptr1,ptr2 

        ! Note: Currently using vx, vy, and vz rather than the adaptive versions

        if (ifld.gt.ldimt1) then  ! For output
          nx1 = lx1u; ny1 = ly1u; nz1 = lz1u
        else
          call getord(ie,ifld)
        end if
     
        if (nx1.ne.lx1) then
          call setmapa(lx1,nx1)
          if (if3d) then
            call specmp(vxe,nx1,u(1,1,1,ie),lx1,ixadpt,ixtadpt,wk1)
            call specmp(vye,ny1,v(1,1,1,ie),ly1,ixadpt,ixtadpt,wk1)
            call specmp(vze,nz1,w(1,1,1,ie),lz1,ixadpt,ixtadpt,wk1)
             
!            CALL MXM (ixadpt,NX1,vx(1,1,1,ie),LX1,wk1,ly1*lz1)
!            DO 5 IZ=1,nz
!              ptr3 = (IZ-1)*NX1*LY1+1
!              CALL MXM (wk1(ptr3),NX1,iytadpt,LY1,wk2,NY1)
! 5          CONTINUE
!            CALL MXM (wk2,lx1*ly1,iztadpt,lz1,vxe,NZ1)
!            ! vy
!            CALL MXM (ixadpt,NX1,vy(1,1,1,ie),LX1,wk1,ly1*lz1)
!            DO 10 IZ=1,nz
!              ptr3 = (IZ-1)*NX1*LY1+1
!              CALL MXM (wk1(ptr3),NX1,iytadpt,LY1,wk2,NY1)
! 10         CONTINUE
!            CALL MXM (wk2,lx1*ly1,iztadpt,lz1,vye,NZ1)
!            ! vz
!            CALL MXM (ixadpt,NX1,vz(1,1,1,ie),LX1,wk1,ly1*lz1)
!            DO 15 IZ=1,nz
!              ptr3 = (IZ-1)*NX1*LY1+1
!              CALL MXM (wk1(ptr3),NX1,iytadpt,LY1,wk2,NY1)
! 15         CONTINUE
!            CALL MXM (wk2,lx1*ly1,iztadpt,lz1,vze,NZ1)
          else
            call specmp(vxe,nx1,u(1,1,1,ie),lx1,ixadpt,ixtadpt,wk1)
            call specmp(vye,ny1,v(1,1,1,ie),ly1,ixadpt,ixtadpt,wk1)
            !! vx
            !CALL MXM (ixadpt,NX1,vx(1,1,1,ie),LX1,wk1,ly1)
            !CALL MXM (wk1,NX1,iytadpt,ly1,vxe,ny1)
            !! vy
            !CALL MXM (ixadpt,NX1,vy(1,1,1,ie),LX1,wk1,ly1)
            !CALL MXM (wk1,NX1,iytadpt,ly1,vye,ny1)
          end if
        else
          nxyz = nx1*ny1*nz1
          call copy(vxe,u(1,1,1,ie),nxyz)
          call copy(vye,v(1,1,1,ie),nxyz)
          if (if3d) call copy(vze,w(1,1,1,ie),nxyz)
        end if
   
      return 
      end 
c-----------------------------------------------------------------------
      subroutine makebdqa
C-----------------------------------------------------------------------
C
C     Add contributions to F from lagged BD terms.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'ADAPT'
C
      COMMON /SCRNSA/ TA (LX1U*LY1U*LZ1U*LELT)
     $ ,              TB (LX1U*LY1U*LZ1U*LELT)
     $ ,              H2 (LX1U*LY1U*LZ1U*LELT)

      integer nel,ntot,ntotall
      integer ptr,ptr1,ptr2 
      real const

C

      NEL   = NELFLD(IFIELD)
      NTOT1 = ntota(ifield)
      CONST = 1./DT
      
      ptr1 = adptr(1,ifield)
      if (ifield.ne.1) ptr = ptr1+ntota(1)
      CALL COPY  (H2,VTRANSA(ptr),NTOT1)
      CALL CMULT (H2,CONST,NTOT1)
C
      CALL COL3  (TB,BM1A(ptr),TAD(ptr1),NTOT1)
      CALL CMULT (TB,BD(2),NTOT1)

      !print *, 'START BDF PART' 
      !print *, 'T', ptr1
      !print '(25F8.3)', tad 
      !print *

      !print *, 'H2'
      !print '(25F8.3)', H2
      !print *
      !print *, 'TB'
      !print '(25F8.3)', TB 
      !print *
C
      ntotall = 0
      do i=2,ifield-1
        ntotall = ntotall + ntota(i)
      end do
      DO 100 ILAG=2,NBD
         !IF (IFGEOM) THEN
         !   CALL COL3 (TA,BM1LAG(1,1,1,1,ILAG-1),
     $   !                 TLAG  (1,1,1,1,ILAG-1,IFIELD-1),NTOT1)
         !ELSE
            ptr2 = ntotall*(lorder-1)+ntot1*(ilag-2)+1
            !print *, 'BDF PTRS', ptr, ptr2
            CALL COL3 (TA,BM1A(ptr),
     $                    TLAGA(ptr2),NTOT1)
            !print *, 'TA', ILAG
            !print '(25F9.6)', TA 
            !print *
            !print *, 'TLAG', ILAG, ptr2
            !print '(25F9.6)', TLAGA 
            !print *
         !ENDIF
         CALL CMULT (TA,BD(ILAG+1),NTOT1)
            !print *, 'TA Again', ILAG
            !print '(25F8.3)', TA 
            !print *
         CALL ADD2  (TB,TA,NTOT1)
 100  CONTINUE
 
      !print *, 'TB Again'
      !print '(25F9.6)', TB 
      !print *
C
      CALL COL2 (TB,H2,NTOT1)
      !print *, 'TB Again 2', ptr1, ntot1
      !print '(25F8.3)', TB 
      !print *
      CALL ADD2 (BQA(ptr1),TB,NTOT1)
      !print *, 'END BDF PART' 
C
      return
      end
c-----------------------------------------------------------------------
      subroutine makeabqa

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'ADAPT'

      COMMON /SCRNZA/ TA (LX1U*LY1U*LZ1U*LELT)
      real AB0, AB1, AB2
      integer NTOT1, NEL, ptr

      AB0   = AB(1)
      AB1   = AB(2)
      AB2   = AB(3)
      NEL   = NELFLD(IFIELD)

      ptr = adptr(1,ifield) 
      ntot1 = ntota(ifield)
      CALL ADD3S2 (TA,VGRADT1A(ptr),
     $                VGRADT2A(ptr),AB1,AB2,NTOT1)
      CALL COPY   (   VGRADT2A(ptr),
     $                VGRADT1A(ptr),NTOT1)
      CALL COPY   (   VGRADT1A(ptr),
     $                     BQA(ptr),NTOT1)
      CALL ADD2S1 (BQA(ptr),TA,AB0,NTOT1)

      return
      end 
c-----------------------------------------------------------------------
      SUBROUTINE FACEISA(CB,S,IEL,IFACE,NX,NY,NZ)
C
C     Assign inflow boundary conditions to face(IE,IFACE)
C     for scalar S.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
      INCLUDE 'TSTEP'     ! ifield    11/19/2010
      INCLUDE 'SOLN'      ! tmask()   11/19/2010
      INCLUDE 'ADAPT'
C
      DIMENSION S(LX1U*LY1U*LZ1U)
      CHARACTER CB*3
c
      common  /nekcb/ cb3
      character*3 cb3
      integer ctr,ptra
      integer ptr,ptr1,ptr2 
      cb3 = cb

      ifld1 = ifield-1


C     Passive scalar term

      ieg = lglel(iel)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)

      IF (CB.EQ.'t  ') THEN
         ptra = adptr(iel,ifield)-1
         call genxyz1(xm1e,ym1e,zm1e,nx,ny,nz,iel)
         NXY = nx*ny
         DO 90 IZ=KZ1,KZ2                           
         DO 90 IY=KY1,KY2                           
         DO 90 IX=KX1,KX2
            ctr = (IZ-1)*NXY+(IY-1)*NX+IX 
            if (tmaska(ptra+ctr).eq.0) then 
               CALL NEKASGN (ctr,1,1,IEL)
               CALL USERBC  (IX,IY,IZ,IFACE,IEG)
               S(ctr) = TEMP
            endif
  90     CONTINUE
        RETURN
C
      ELSEIF (CB.EQ.'ms ' .OR. CB.EQ.'msi') THEN
C
         DO 200 IZ=KZ1,KZ2
         DO 200 IY=KY1,KY2
         DO 200 IX=KX1,KX2
            ctr = (IZ-1)*NXY+(IY-1)*NX+IX 
            CALL NEKASGN (ctr,1,1,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            S(ctr) = SIGMA
  200    CONTINUE
C
!      ELSEIF (CB.EQ.'kd ') THEN
!C
!         DO 300 IZ=KZ1,KZ2
!         DO 300 IY=KY1,KY2
!         DO 300 IX=KX1,KX2
!            CALL NEKASGN (IX,IY,IZ,IEL)
!            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
!            S(IX,IY,IZ) = TURBK
!  300    CONTINUE
!C
!      ELSEIF (CB.EQ.'ed ') THEN
!C
!         DO 400 IZ=KZ1,KZ2
!         DO 400 IY=KY1,KY2
!         DO 400 IX=KX1,KX2
!            CALL NEKASGN (IX,IY,IZ,IEL)
!            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
!            S(IX,IY,IZ) = TURBE
!  400    CONTINUE
!C
      ENDIF
C
      RETURN
      END
c----------------------------------------------------------------------
      subroutine sv3d_modal(glo_num,ngv,nx,ny,nz,nel,vertex,ifcenter)

c     setup unique ids for dssum  

c     note:
c     total number of unique vertices, edges and faces has to be smaller 
c     than 2**31 (integer-4 limit).
c     if nelgt < 2**31/12 we're ok for sure (independent of N)! 

      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'

      integer*8 glo_num(1),ngv
      integer vertex(0:1,0:1,0:1,nel),nx(nel),ny(nel),nz(nel)
      integer nx27(27,nel),ny27(27,nel),nz27(27,nel)
      logical ifcenter

      integer  edge(0:1,0:1,0:1,3,lelt),enum(12,lelt),fnum(6,lelt)
      common  /scrmg/ edge,enum,fnum

      parameter (nsafe=8)  ! OFTEN, nsafe=2 suffices
      integer etuple(4,12*lelt*nsafe),ftuple(5,6,lelt*nsafe)
      integer ind(12*lelt*nsafe)
      common  /scrns/ ind,etuple
      equivalence  (etuple,ftuple)

      integer gvf(4),facet(4),aa(3),key(3),e
      logical ifij
      
      integer*8 igv,ig0
      integer*8 ngvv,ngve,ngvs,ngvi,ngvm
      integer*8 n_on_edge,n_on_face,n_in_interior
      integer*8 i8glmax

      key(1)=1
      key(2)=2
      key(3)=3

      call get_nx27(nx27,ny27,nz27,nx,ny,nz,nel)

c     Assign hypercube ordering of vertices
c     -------------------------------------

c     Count number of unique vertices
      nlv  = 2**ndim
      ngvv = iglmax(vertex,nlv*nel)

      nxu = iglmax(nx,nel)
      nyu = iglmax(ny,nel)
      nzu = iglmax(nz,nel)
      nu  = max(nxu,nyu)
      nu  = max(nu,nzu)

      ioff = 0
      do e=1,nel
         do k=0,1
         do j=0,1
         do i=0,1
c           Local to global node number (vertex)
            il =1+(nx(e)-1)*i+nx(e)*(ny(e)-1)*j+nx(e)*ny(e)*(nz(e)-1)*k
            ile=il + ioff
            glo_num(ile)   = vertex(i,j,k,e)
         enddo
         enddo
         enddo
         ioff = ioff + nx(e)*ny(e)*nz(e)
      enddo
      ngv  = ngvv

      if (nu.eq.2) return

c     Assign global vertex numbers to SEM nodes on each edge
c     ------------------------------------------------------

      do e=1,nel  ! Assign edge labels by bounding vertices.  
         do k=0,1
         do j=0,1
         do i=0,1
            edge(i,j,k,1,e) = vertex(i,j,k,e)  ! r-edge
            edge(j,i,k,2,e) = vertex(i,j,k,e)  ! s-edge
            edge(k,i,j,3,e) = vertex(i,j,k,e)  ! t-edge
         enddo
         enddo
         enddo
      enddo

      do i=0,12*nel-1 ! Sort edges by bounding vertices.
         if (edge(0,i,0,1,1).gt.edge(1,i,0,1,1)) then
            kswap = edge(0,i,0,1,1)
            edge(0,i,0,1,1) = edge(1,i,0,1,1)
            edge(1,i,0,1,1) = kswap
         endif
         etuple(3,i+1) = edge(0,i,0,1,1)
         etuple(4,i+1) = edge(1,i,0,1,1)
      enddo

      m    = 4              ! Assign a number (rank) to each unique edge
      n    = 12*nel
      nmax = 12*lelt*nsafe  ! nsafe for crystal router factor of safety
      call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
      do i=1,12*nel
         enum(i,1) = etuple(3,i)
      enddo
      n_unique_edges = iglmax(enum,12*nel)

      !print *, 'GLO_NUM 1'
      !print '(16I8)', glo_num
      !print *

      n_on_edge = nu-2
      ngve      = n_unique_edges*n_on_edge
      ie = 0
      do e=1,nel
         iedg_loc = 0

         do k=0,1  !  Edges 1-4
         do j=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = nx(e)*(ny(e)-1)*j + nx(e)*ny(e)*(nz(e)-1)*k
            i0e = i0 + ie
            i27  = 2*(1 + j*3 + k*9)
            do i=2,nx27(i27,e)-1
               glo_num(i0e+i) = igv + i-1
            enddo
            iedg_loc = iedg_loc + 1
         enddo
         enddo

         do k=0,1  ! Edges 5-8
         do i=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = 1+(nx(e)-1)*i + nx(e)*ny(e)*(nz(e)-1)*k
            i0e = i0 + ie
            i27  = 2*(2 + i + k*9)
            do j=2,ny27(i27,e)-1
               glo_num(i0e+(j-1)*nx(e)) = igv + j-1
            enddo
            iedg_loc = iedg_loc + 1
         enddo
         enddo
c
c        Edges 9-12
         do j=0,1
         do i=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = 1 + (nx(e)-1)*i + nx(e)*(ny(e)-1)*j
            i0e = i0 + ie
            i27  = 2*(5 + i + j*3)
            do k=2,nz27(i27,e)-1
               glo_num(i0e+(k-1)*nx(e)*ny(e)) = igv + k-1
            enddo
            iedg_loc = iedg_loc + 1
         enddo
         enddo
         ie = ie + nx(e)*ny(e)*nz(e)
      enddo
      ngv   = ngv + ngve
c
c     Assign global node numbers on the interior of each face
c     ------------------------------------------------------ 
c
c     Assign faces by 3-tuples 
c
c     (The following variables all take the symmetric 
c     notation of IFACE as arguments:)
c
c     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
c                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.
c
c                        3+-----+4    ^ Y
c                        /  2  /|     |
c     Edge 1 extends    /     / |     |
c       from vertex   7+-----+8 +2    +----> X
c       1 to 2.        |  4  | /     /
c                      |     |/     /
c                     5+-----+6    Z
c                         3
c
      nfaces=ndim*2
      ncrnr =2**(ndim-1)
      do e=1,nel
         do ifac=1,nfaces
            do icrn=1,ncrnr
               i                  = icface(icrn,ifac)-1
               facet(icrn)        = vertex(i,0,0,e)
            enddo
            call isort(facet,ind,ncrnr)
            call icopy(ftuple(3,ifac,e),facet,ncrnr-1)
         enddo
      enddo

c     Assign a number (rank) to each unique face
      m    = 5
      n    = 6*nel
      nmax = 6*lelt*nsafe  ! nsafe for crystal router factor of safety
      call gbtuple_rank(ftuple,m,n,nmax,cr_h,nid,np,ind)
      do i=1,6*nel
         fnum(i,1) = ftuple(3,i,1)
      enddo
      n_unique_faces = iglmax(fnum,6*nel)

      !print *, 'GLO_NUM 2'
      !print '(16I8)', glo_num
      !print *

      ie = 0
      do e=1,nel
       call dsset (nx(e),ny(e),nz(e))
       do iface=1,nfaces
         i0 = skpdat(1,iface)
         i1 = skpdat(2,iface)
         is = skpdat(3,iface)
         j0 = skpdat(4,iface)
         j1 = skpdat(5,iface)
         js = skpdat(6,iface)

c        On each face, count from minimum global vertex number,
c        towards smallest adjacent vertex number.  e.g., suppose
c        the face is defined by the following global vertex numbers:


c                    11+--------+81
c                      |c      d|
c                      |        |
c                      |        |
c                      |a      b|
c                    15+--------+62
                          
c        We would count from c-->a, then towards d.

         if (iface==1.or.iface==3) then
            nr = nx(e); ns = ny(e)
         elseif (iface==2.or.iface==4) then
            nr = nx(e); ns = ny(e)
         else
            nr = nx(e); ns = ny(e)
         end if

         gvf(1) = glo_num(i0+nr*(j0-1)+ie)
         gvf(2) = glo_num(i1+nr*(j0-1)+ie)
         gvf(3) = glo_num(i0+nr*(j1-1)+ie)
         gvf(4) = glo_num(i1+nr*(j1-1)+ie)

         call irank(gvf,ind,4)

c        ind(1) tells which element of gvf() is smallest.

         ifij = .false.
         if (ind(1).eq.1) then
            idir =  1
            jdir =  1
            if (gvf(2).lt.gvf(3)) ifij = .true.
         elseif (ind(1).eq.2) then
            idir = -1
            jdir =  1
            if (gvf(1).lt.gvf(4)) ifij = .true.
         elseif (ind(1).eq.3) then
            idir =  1
            jdir = -1
            if (gvf(4).lt.gvf(1)) ifij = .true.
         elseif (ind(1).eq.4) then
            idir = -1
            jdir = -1
            if (gvf(3).lt.gvf(2)) ifij = .true.
         endif
c
         if (idir.lt.0) then
            it=i0
            i0=i1
            i1=it
            is=-is
         endif
c
         if (jdir.lt.0) then
            jt=j0
            j0=j1
            j1=jt
            js=-js
         endif
c
         nxx = nr*ns
         n_on_face = (nr-2)*(ns-2) 
         ngvs  = n_unique_faces*n_on_face
         ig0 = ngv + n_on_face*(fnum(iface,e)-1)
         if (ifij) then
            k=0
            l=0
            do j=j0,j1,js
            do i=i0,i1,is
               !print *, e,iface,i,j,i+nr*(j-1)
               k=k+1
c              this is a serious kludge to stay on the face interior
               if (k.gt.nr.and.k.lt.nxx-nr .and.
     $            mod(k,nr).ne.1.and.mod(k,nr).ne.0) then
c                 interior
                  l = l+1
                  glo_num(i+nr*(j-1)+ie) = l + ig0
               endif
            enddo
            enddo
         else
            k=0
            l=0
            do i=i0,i1,is
            do j=j0,j1,js
               k=k+1
c              this is a serious kludge to stay on the face interior
               if (k.gt.nr.and.k.lt.nxx-nr .and.
     $            mod(k,nr).ne.1.and.mod(k,nr).ne.0) then
c                 interior
                  l = l+1
                  glo_num(i+nr*(j-1)+ie) = l + ig0
               endif
            enddo
            enddo
         endif
         ! Change node numbering if backwards
         ! rows
!         print *, ie+i0+2*is, ie+i0+is
!         if (glo_num(ie+i0+2*is)<glo_num(ie+i0+is)) then  ! Need to somehow account for the zero location
!           call FaceReverseRows(glo_num(ie+1),nr,ns,iface)
!         end if 
!         ! columns
!         print *, ie+i0+is+nr*(j0+2*js-1),ie+i0+is+nr*(j0+js-1)
!
!         if (glo_num(ie+i0+is+nr*(j0+2*js-1))
!     $       <glo_num(ie+i0+is+nr*(j0+js-1))) then
!           call FaceReverseCols(glo_num(ie+1),nr,ns,iface)
!         end if 
         
       enddo
       ie = ie + nx(e)*ny(e)*nz(e)
 
      enddo
      ngv   = ngv + ngvs
      
      !print *, 'GLO_NUM 2'
      !print '(16I8)', glo_num
      !print *
c
c     Finally, number interiors (only ifcenter=.true.)
c     -------------------------------------------------
c
      n_in_interior = (nx(e)-2)*(ny(e)-2)*(nz(e)-2)
      ngvi = n_in_interior*nelgt
      if (ifcenter) then
         do e=1,nel
            ig0 = ngv + n_in_interior*(lglel(e)-1)
            l = 0
            do k=2,nz(e)-1
            do j=2,ny(e)-1
            do i=2,nx(e)-1
               l = l+1
               glo_num(i+nx(e)*(j-1)
     $                 +nx(e)*ny(e)*(k-1)+nxyz*(e-1)) = ig0+l
            enddo
            enddo
            enddo
         enddo
         ngv = ngv + ngvi
      else
         do e=1,nel
            l = 0
            do k=2,nz(e)-1
            do j=2,ny(e)-1
            do i=2,nx(e)-1
               l = l+1
               glo_num(i+nx(e)*(j-1)+nx(e)*ny(e)*(k-1)+nxyz*(e-1)) = 0
            enddo
            enddo
            enddo
         enddo
      endif
c
c     Quick check on maximum #dofs:
!      m    = nxyz*nelt
!      ngvm = i8glmax(glo_num,m)
!      ngvv = ngvv + ngve + ngvs  ! number of unique ids w/o interior 
!      ngvi = ngvi + ngvv         ! total number of unique ids 
!      if (nid.eq.0) write(6,1) nx,ngvv,ngvi,ngv,ngvm
!    1 format('   setvert3d:',i4,4i12)
c
      return
      end
c------------------------------------------------------------------------
      subroutine get_nx27(nx27,ny27,nz27,nx,ny,nz,nel)
 
      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
 
      real nx27r(27,nel),ny27r(27,nel),nz27r(27,nel),nxr,nyr,nzr
      integer nx(nel),ny(nel),nz(nel),
     $        nx27(27,nel),ny27(27,nel),nz27(27,nel)
      integer ifld
      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex
      integer*8 gnum27(27*nel),ngv

      ! Fill 3x3 cube with order
      do ie=1,nel
        nxr = real(nx(ie))
        nyr = real(ny(ie))
        nzr = real(nz(ie))
        do i=1,27
           nx27r(i,ie) = nxr
           ny27r(i,ie) = nyr
           nz27r(i,ie) = nzr
        end do 
      end do

      ! Get minimum
      call setupds(gsh_fld(ldimt1+1),3,3,3,nelt,nelgt,vertex,gnum27)

      ifld = ifield
      ifield = ldimt1+1
      call dsop(nx27r,'m  ',3,3,3)
      call dsop(ny27r,'m  ',3,3,3)
      call dsop(nz27r,'m  ',3,3,3)
      ifield = ifld

      ! convert back to integer
      do ie=1,nel
        do i=1,27
           nx27(i,ie) = int(nx27r(i,ie))
           ny27(i,ie) = int(ny27r(i,ie))
           nz27(i,ie) = int(nz27r(i,ie))
        end do 
      end do

      return 
      end
c-----------------------------------------------------------------------
      subroutine FaceReverseRows(gnum,nr,iface)

      include 'SIZE' 
      include 'TOPOL' 

      integer gnum(1),nr,iface,newrow(nr),ctr

      i0 = skpdat(1,iface)
      i1 = skpdat(2,iface)
      is = skpdat(3,iface)
      j0 = skpdat(4,iface)
      j1 = skpdat(5,iface)
      js = skpdat(6,iface)

      do j=j0,j1,js
       ctr = 0
       do i=i0,i1,is
         ctr = ctr+1
         newrow(ctr) = gnum(i+nr*(j-1)-1)
       end do
       ctr = 0
       do i=i1,i0,-is
         ctr = ctr+1
         gnum(i+nr*(j-1)-1) = newrow(ctr)
       end do
      end do

      return
      end 
c-----------------------------------------------------------------------
      subroutine FaceReverseCols(gnum,nr,iface)
      
      include 'SIZE' 
      include 'TOPOL' 

      integer gnum(1),nr,iface,newcol(nr),ctr

      i0 = skpdat(1,iface)
      i1 = skpdat(2,iface)
      is = skpdat(3,iface)
      j0 = skpdat(4,iface)
      j1 = skpdat(5,iface)
      js = skpdat(6,iface)

      do i=i0,i1,is
       ctr = 0
       do j=j0,j1,js
         ctr = ctr+1
         newcol(ctr) = gnum(i+nr*(j-1)-1)
       end do
       ctr = 0
       do j=j1,j0,-js
         ctr = ctr+1
         gnum(i+nr*(j-1)-1) = newcol(ctr)
       end do
      end do

      return
      end 
c-----------------------------------------------------------------------
C      subroutine dssum_modal(u,nx,ny,nz)
C
C      include 'SIZE'
C      include 'ADAPT'
C      include 'TSTEP'
C
C      real u(1)
C      integer nx,ny,nz      
C      real wk1(lx1u),wk2(lx1u*ly1u),wk3(lx1u*ly1u*lz1u),
C     $     uh(lx1u*ly1u*lz1u*lelt),Lj(lx1u*lx1u),Ljt(lx1u*lx1u)
C
C      if (minord(ifield).ne.maxord(ifield)) then
C        ! Convert to modal space
C        nel = nelfld(ifield)
C        do ie=1,nel
C          ptr = adptr(ie,ifield)
C          call zwgll(wk1,wk2,nx1)  ! wk1 holds the nodes
C          call build_legend_transform(Lj,Ljt,wk1,nx1)  ! NOTE: Need to use different basis ordering
C          call tensr3(uh(ptr),nx,u(ptr),nx,Lj,Ljt,Ljt,wk3)    !  Go to Legendre space
C        end do      
C
C        ! dssum
C        call dssum(uh,nx,ny,nz)
C
C        ! Convert back to nodal space
C        do ie=1,nel
C          ptr = adptr(ie,ifield)
C          call zwgll(wk1,wk2,nx1)  ! wk1 holds the nodes
C          call build_legend(Lj,Ljt,wk1,nx1)
C          call tensr3(uh(ptr),nx,u(ptr),nx,Lj,Ljt,Ljt,wk3)    !  Go to nodal space
C        end do      
C      else
C        call dssum(u,nx,ny,nz)   
C      end if   
C
C      return
C      end  
c-----------------------------------------------------------------------
      subroutine build_legend(Lj,Ljt,zpts,nx)
c
      real Lj(nx*nx),Ljt(nx*nx),zpts(nx)
c
      parameter (lm=90)
      integer   indr(lm),indc(lm),ipiv(lm),nx

      j = 1
      n = nx-1
      do i=1,nx
         z = zpts(i)
         call legendre_poly(Lj(j),z,n)  ! Return Lk(z), k=0,...,n
         j = j+nx
      enddo
      call transpose1(Lj,nx)
      call transpose (Ljt,nx,Lj,nx)
      
      return
      end
c-----------------------------------------------------------------------
C      subroutine UpdateXYZ
C     
C      include 'SIZE' 
C      include 'ADAPT'
C      include 'GEOM'
C     
C      integer ifld,ie,ntot,nxyz
C 
C      ntot = ivlsum(ntota,ldimt1)
C      call copy(soltmp,bm1a,ntot) 
C      do ifld=1,ldimt1
C        do ie=1,nelt
C          ptr1 = adptr(ie,ifld) 
C          ptr2 = newptr(ie,ifld)
C          if (ifld.ne.1) then
C             ptr1 = ptr1+ntota(1)
C             ptr2 = ptr2+newntota(1)
C          endif
C          if (adaptflag(ie,ifld)) then
C            call getord(ie,ifld) 
C            nx1 = nx1+adinc; ny1 = ny1+adinc;
C            if (nz1.ne.1) nz1 = nz1+adinc
C            nxyz = nx1*ny1*nz1
C            call genxyz1(xm1a(ptr2),ym1a(ptr2),zm1a(ptr2),
C     $                   nx1,ny1,nz1,ie)
C          else
C            nxyz = nx1*ny1*nz1
C            call copy(xm1a(ptr2),soltmp(ptr1),nxyz)
C            call copy(ym1a(ptr2),soltmp(ptr1),nxyz)
C            call copy(zm1a(ptr2),soltmp(ptr1),nxyz)
C          end if
C        end do
C      end do 
C
C      end
c------------------------------------------------------------------
Cc------------------------------------------------------------------
C      subroutine setupdsa(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
C
C      include 'SIZE'
C      include 'INPUT'
C      include 'PARALLEL'
C      include 'NONCON'
C
C      integer gs_handle
C      integer vertex(1)
C      integer*8 glo_num(1),ngv
C
C      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
C
C
Cc     Global-to-local mapping for gs
C      call set_verta(glo_num,ngv,nx,nel,vertex,.false.)
C
Cc     Initialize gather-scatter code
C      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)
C
Cc     call gs_chkr(glo_num)
C
C      return
C      end
Cc------------------------------------------------------------------
C      subroutine set_verta(glo_num,ngv,nx,nel,vertex,ifcenter)
Cc
Cc     Given global array, vertex, pointing to hex vertices, set up
Cc     a new array of global pointers for an nx^ndim set of elements.
Cc
C      include 'SIZE'
C      include 'INPUT'
Cc
C      integer*8 glo_num(1),ngv
C      integer vertex(1),nx
C      logical ifcenter
C
C      if (if3d) then
C         call setvert3da(glo_num,ngv,nel,vertex,ifcenter)
C      else
C         call setvert2da(glo_num,ngv,nel,vertex,ifcenter)
C      endif
C
C      ! NOT ADAPTIVE: not sure what this does below so it
C      ! is currently not part of the adaptive step
C      
Cc     Check for single-element periodicity 'p' bc
C      !nz = 1
C      !if (if3d) nz = nx
C      !call check_p_bc(glo_num,nx,nx,nz,nel)
C
C      !if(nid.eq.0) write(6,*) 'call usrsetvert'
C      !call usrsetvert(glo_num,nel,nx,nx,nx)
C      !if(nid.eq.0) write(6,'(A,/)') ' done :: usrsetvert'
C
C      return
C      end
Cc------------------------------------------------------------------
C      subroutine setvert2da(glo_num,ngv,nel,vertex,ifcenter)
Cc
Cc     setup unique ids for dssum  
Cc
C      include 'SIZE'
C      include 'CTIMER'
C      include 'PARALLEL'
C      include 'TOPOL'
C      include 'GEOM'
C      include 'ADAPT'
C
C      integer*8 glo_num(1),ngv
C      integer vertex(0:1,0:1,1)
C      logical ifcenter
C
C      integer  edge(0:1,0:1,2,lelt),enum(4,lelt)
C      common  /scrmg/ edge,enum
C
C      parameter (nsafe=8)  ! OFTEN, nsafe=2 suffices
C      integer etuple(4,4*lelt*nsafe),ind(4*lelt*nsafe)
C      common  /scrns/ ind,etuple
C
C      integer gvf(4),aa(3),key(3),e,eg
C      logical ifij
C
C      integer*8 igv,ig0
C      integer*8 ngvv,ngve,ngvs,ngvi,ngvm
C      integer*8 n_on_edge,n_on_face,n_in_interior
C      integer*8 i8glmax
Cc
Cc
C      key(1)=1
C      key(2)=2
C      key(3)=3
Cc
Cc     Count number of unique vertices
C      nlv  = 2**ndim
C      ngvv = iglmax(vertex,nlv*nel)
C      ngv  = ngvv
Cc
Cc     Assign hypercube ordering of vertices.
C      do e=1,nel
C         call getord(e,ifield)
C         do j=0,1
C         do i=0,1
Cc           Local to global node number (vertex)
C            il  = (nx1-1)*i + nx1*(ny1-1)*j
C            ile = adptr1(e,ifield)+il
C            glo_num(ile)   = vertex(i,j,e)
C         enddo
C         enddo
C      enddo
C      if (nx.eq.2) return
Cc
Cc     Assign edge labels by bounding vertices.  
C      do e=1,nel
C         do j=0,1
C         do i=0,1
C            edge(i,j,1,e) = vertex(i,j,e)  ! r-edge
C            edge(j,i,2,e) = vertex(i,j,e)  ! s-edge
C         enddo
C         enddo
C      enddo
C
Cc     Sort edges by bounding vertices.
C      do i=0,4*nel-1
C         if (edge(0,i,1,1).gt.edge(1,i,1,1)) then
C            kswap = edge(0,i,1,1)
C            edge(0,i,1,1) = edge(1,i,1,1)
C            edge(1,i,1,1) = kswap
C         endif
C         etuple(3,i+1) = edge(0,i,1,1)
C         etuple(4,i+1) = edge(1,i,1,1)
C      enddo
C
Cc     Assign a number (rank) to each unique edge
C      m    = 4
C      n    = 4*nel
C      nmax = 4*lelt*nsafe  ! nsafe for crystal router factor of safety
C
C      call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
C      do i=1,4*nel
C         enum(i,1) = etuple(3,i)
C      enddo
C      n_unique_edges = iglmax(enum,4*nel)
C
Cc= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
Cc     Assign global vertex numbers to SEM nodes on each edge
C
C      call dsop(adords(1,1,1,1,ifield),"min",3,3,3)
C
C      do e=1,nel
C         iedg_loc = 0
C         call getord(e,ifield)
C         n_on_edge = nx1-2
Cc        Edges 1-2
C         do j=0,1
C            !!!!!! NEED TO FIX IGV
C            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1) ! number of previous global vertices 
C            i0  = nx1*(ny1-1)*j                            ! local position  
C            i0e = i0 + adptr1(e,ifield)-1                        ! global position
C            if (glo_num(i0e+1).lt.glo_num(i0e+nx1)) then 
C               do i=2,nx-1                                   ! std forward case
C                  glo_num(i0e+i) = igv + i-1
C               enddo
C            else
C               do i=2,nx-1                                   ! backward case
C                  glo_num(i0e+i) = igv + 1 + n_on_edge-(i-1)  ! Need to put in the appropriate n_on_edge
C               enddo
C            endif
C            iedg_loc = iedg_loc + 1
C         enddo
Cc
Cc        Edges 3-4
C         do i=0,1
C            n_on_edge = ny1-2
C            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
C            i0  = 1+(nx1-1)*i
C            i0e = i0 + adptr1(e,ifield)
C            if (glo_num(i0e).lt.glo_num(i0e+nx1*(ny1-1))) then
C               do j=2,ny1-1                                   ! std forward case
C                  glo_num(i0e+(j-1)*nx1) = igv + j-1
C               enddo
C            else
C               do j=2,ny1-1                                   ! backward case
C                  glo_num(i0e+(j-1)*nx1) = igv + 1 + n_on_edge-(j-1)
C               enddo
C            endif
C            iedg_loc = iedg_loc + 1
C         enddo
C      enddo
C 
C      !ngve = n_unique_edges*n_on_edge 
C      !ngv  = ngv + ngve    
Cc
Cc     Finally,  number interiors  
Cc
C      !n_in_interior = (nx1-2)*(ny1-2)
C      !ngvi          = n_in_interior*nelgt
C      if (ifcenter) then
C        ! do e=1,nel
C        !    call getord(e,ifield)
C        !    n_in_interior = (nx1-2)*(ny1-2)
C        !    ig0 = ngv + n_in_interior*(lglel(e)-1)
C        !    l = 0
C        !    ptr = adptr1(e,ifield)-1
C        !    do j=2,ny1-1
C        !    do i=2,nx1-1
C        !       l = l+1
C        !       glo_num(i+nx1*(j-1)+ptr) = ig0+l
C        !    enddo
C        !    enddo
C        ! enddo
C        ! ngv = ngv + ngvi
C      else
C         do e=1,nel
C            call getord(e,ifield)
C            n_in_interior = (nx1-2)*(ny1-2)
C            l = 0
C            ptr = adptr1(e,ifield)-1
C            do j=2,ny1-1
C            do i=2,nx1-1
C               l = l+1
C               glo_num(i+nx1*(j-1)+ptr) = 0
C            enddo
C            enddo
C         enddo
C      endif
C
C      return
C      end
c------------------------------------------------------------------
!      subroutine errest(uerr,terr)
! 
!      real nvals
!
!      ! Section 7.4.4 in Deville, Fischer, and Mund covers the error estimator coded here
!
!      !do i=1,ldimt
!      do e = 1,nelt
!        ! Get coefficients for Legendre basis
!        call chofbasis(coef)
!
!        ! Do a least squares fit of a(n) = C*exp(-alpha*n) for the last 4 coefficients
!        if (nx1tva>3) then
!          num = 4
!        else
!          num = nx1
!        end if
!        nvals = real((/nx1-num+1:nx1:1/))
!        call lsqexp(nvals,coef(nx1-num+1),num,terr(e,1),C)  ! terr = alpha
!        
!        ! if alpha < 1, increase order or just return alpha as terr and change inequatily in calling routine
!
!        !! consider ways of decreasing order as well (for alpha >> 1)
!      end do
!      !end do
!
!      end
!c-----------------------------------------------------------------------
!      subroutine lsqexp(x,y,n,alpha,C)
! 
!      ! least squares fit for y = C*exp(alpha*x) where n the length of the data
! 
!      integer i, n
!      real y(*), x(*), alpha, C
!      real sumy,sumxy,sumx2y,sumylny,sumxylny,lny,xi,yi
!      
!      sumx = 0.0; sumx2 = 0.0; sumlny = 0.0; sumxlny = 0.0
!      do i=1,n
!        xi = x(i)
!        yi = y(i)
!        lny = alog(yi)
!        sumy = sumy+yi
!        sumxy = sumxy+xi*yi
!        sumx2y = sumx2y+xi*xi*yi
!        sumylny = sumylny+yi*lny 
!        sumxylny = sumxylny+xi*yi*lny
!      end do
! 
!      C = 1.0/(sumy*sumx2y-sumxy*sumxy) 
!      alpha = (sumy*sumxylny-sumxy*sumylny)*C
!      C = (sumx2y*sumylny-sumxy*sumxylny)*C
!      C = exp(C)
!
!      end 
c-----------------------------------------------------------------------
