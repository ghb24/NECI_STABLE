!TODO: Make sure everything has an ierr argument when allocating, and that everything is deallocated at end
!TODO: Add symmetry information by using singles excitation generators (store excitations)? (Not urgent) 
!Module to non-iteratively calculate the RPA energy under the quasi-boson approximation
module RPA_Mod
    use SystemData, only: nel, nBasis, Arr, Brr, G1, tReltvy
    use sltcnd_mod, only: sltcnd_2
    use constants, only: dp, int64, n_int
    use Determinants, only: get_helement, fDet
    use SymExcit3, only: GenExcitations3
    use SymExcit4, only: GenExcitations4, ExcitGenSessionType
    use Determinants, only: GetH0Element3
    use bit_reps, only: NIfTot
    use DetBitops, only: EncodeBitDet
    use Integrals_neci, only: get_umat_el
    use UMatCache, only: GTID
    use util_mod, only: get_free_unit

    implicit none

    !input options
    logical :: tStabilityAnalysis=.true. 
    logical :: tDirectRPA

    !global variables
    integer :: ov_space
    integer :: virt_start

    !global arrays
    real(dp), allocatable :: A_mat(:,:)
    real(dp), allocatable :: B_mat(:,:)
    real(dp), allocatable :: OccNumbers(:)
    
    contains

    subroutine RunRPA_QBA(Weight,Energy)
        implicit none
        integer :: ierr,i,j,m,n,ex(2,2),ex2(2,2),mi_ind,nj_ind
        integer :: StabilitySize,lWork,info,i_p,m_p,v,mp_ip_ind,ic
        integer :: nJ(NEl),exflag,mu,id(2,2),a,i_ind,iunit,ia_ind
        integer(n_int) :: iLutHF(0:NIfTot)
        real(dp), intent(out) :: Weight,Energy
        real(dp) :: Energy_stab,Temp_real,norm,Energy2,H0tmp,Fii
        real(dp) :: X_norm,Y_norm
        HElement_t(dp) :: HDiagTemp,hel,hel1,hel2
        logical :: tAllExcitsFound,tParity
        real(dp), allocatable :: Stability(:,:),temp2(:,:),W2(:)
        real(dp), allocatable :: W(:),Work(:),S_half(:,:),temp(:,:)
        real(dp), allocatable :: X_stab(:,:),Y_stab(:,:),StabilityCopy(:,:)
        real(dp), allocatable :: X_chol(:,:),Y_chol(:,:)
        real(dp), allocatable :: AminB(:,:),AplusB(:,:),temp3(:,:)
        character(len=*), parameter :: t_r="RunRPA_QBA"

        type(ExcitGenSessionType) :: session

        write(6,"(A)")
        write(6,"(A)") "**************************************"
        if(tDirectRPA) then
        write(6,"(A)") "*   Entering DIRECT RPA calculation  *"
        else
        write(6,"(A)") "*   Entering FULL RPA calculation    *"
        endif
        write(6,"(A)") "**************************************"
        write(6,"(A)") 
        call neci_flush(6)

        HDiagTemp = get_helement(fDet, fDet, 0)
        Energy = real(HDiagTemp,dp)
        write(6,"(A,G25.10)") "Reference energy is: ",Energy

        !Quickly find correlation energy from MP2 for comparison.
        Temp_real = 0.0_dp
        tAllExcitsFound=.false.
        exflag=3
        ex(:,:)=0
        call EncodeBitDet(FDet,iLutHF)
        HDiagTemp=GetH0Element3(FDet)
        Fii=real(HDiagTemp,dp)

        do while(.true.)
            if (tReltvy) then 
                call GenExcitations4(session, FDet, nJ, exFlag, Ex, tParity, tAllExcitsFound, .false.)
            else
                call GenExcitations3(FDet,iLutHF,nJ,exflag,Ex,tParity,tAllExcitsFound,.false.)
            endif
            if(tAllExcitsFound) exit !All excits found
            if(Ex(1,2).eq.0) then
                ic=1
            else
                ic=2
            endif
            hel=get_helement(FDet,nJ,ic,Ex,tParity)
            H0tmp=getH0Element3(nJ)
            H0tmp=Fii-H0tmp
            Temp_real=Temp_real+(hel**2)/H0tmp
        enddo
        write(6,"(A,G25.10)") "For comparison, MP2 correlation energy is: ",Temp_real

        Weight=(0.0_dp)
        Energy=(0.0_dp)

        ov_space=(nBasis-NEl)*NEl
        virt_start=NEl+1
        write(6,"(A,I8)") "1p-1h space is: ",ov_space

        allocate(A_mat(ov_space,ov_space),stat=ierr)
        allocate(B_mat(ov_space,ov_space),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc Err")
        A_mat(:,:)=0.0_dp
        B_mat(:,:)=0.0_dp

        !First, construct A and B, which correspond to:
        !  <HF| [a*_i a_m [H, a*_n a_j]] |HF> = A
        ! -<HF| [a*_i a_m [H, a*_j a_n]] |HF> = B
        do j=1,nel
            ex(1,2)=Brr(j)   !Second index in integral
            ex2(2,2)=Brr(j)
            do n=virt_start,nBasis
                nj_ind = ov_space_ind(n,j)
                ex(2,2)=Brr(n)   !fourth index in integral
                ex2(1,2)=Brr(n)
                do i=1,nel
                    ex(2,1)=Brr(i)   !Third index in integral
                    ex2(2,1)=Brr(i)  
                    do m=virt_start,nBasis
                        mi_ind = ov_space_ind(m,i)
                        ex(1,1)=Brr(m)   !First index in integral
                        ex2(1,1)=Brr(m) 
 
                        if(tDirectRPA) then
                            !No exchange interactions
                            ! Obtain spatial rather than spin indices if required
                            id = gtID(ex)
                            ! Only non-zero contributions if Ms preserved in each term (consider
                            ! physical notation).
                            if ( (G1(ex(1,1))%Ms == G1(ex(2,1))%Ms) .and.                   &
                                 (G1(ex(1,2))%Ms == G1(ex(2,2))%Ms) ) then
                                hel1 = get_umat_el (id(1,1), id(1,2), id(2,1), &
                                                    id(2,2))
                            else
                                hel1 = (0)
                            endif
                            A_mat(mi_ind,nj_ind) = real(hel1,dp)
                            !Now for B matrix
                            id = gtID(ex2)
                            ! Only non-zero contributions if Ms preserved in each term (consider
                            ! physical notation).
                            if ( (G1(ex2(1,1))%Ms == G1(ex2(2,1))%Ms) .and. &
                                 (G1(ex2(1,2))%Ms == G1(ex2(2,2))%Ms) ) then
                                hel1 = get_umat_el (id(1,1), id(1,2), id(2,1), &
                                                    id(2,2))
                            else
                                hel1 = (0)
                            endif
                            B_mat(mi_ind,nj_ind) = real(hel1,dp)
                        else
                            !Full antisymmetrized integrals
                            HEl1 = sltcnd_2(ex,.false.)
                            HEl2 = sltcnd_2(ex2,.false.)
                            A_mat(mi_ind,nj_ind) = real(HEl1,dp)           
                            B_mat(mi_ind,nj_ind) = real(HEl2,dp)           
                        endif

                    enddo
                enddo
            enddo
        enddo

        !Now add the diagonal part to A
        do i=1,nel
            do m=virt_start,nBasis
                mi_ind = ov_space_ind(m,i)

                A_mat(mi_ind,mi_ind) = A_mat(mi_ind,mi_ind) + (Arr(Brr(m),2)-Arr(Brr(i),2))
!                write(6,*) Arr(Brr(m),2)-Arr(Brr(i),2)
            enddo
        enddo

        !Check that A is hermitian and B is symmetric
        do i=1,ov_space
            do j=1,ov_space
                if(abs(B_mat(i,j)-B_mat(j,i)).gt.1.0e-7_dp) then
                    write(6,*) i,j,B_mat(i,j), B_mat(j,i),abs(B_mat(i,j)-B_mat(j,i))
                    call stop_all(t_r,"B not symmetric")
                endif
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7_dp) then
                    write(6,*) i,j,A_mat(i,j), A_mat(j,i),abs(A_mat(i,j)-A_mat(j,i))
                    
                    call stop_all(t_r,"A not hermitian")
                endif
            enddo
        enddo

        if(tStabilityAnalysis) then
            !Calculate stability analysis of HF solution, and from this, calculate RPA amplitudes
            !This should give identical results to other method, will be quite a bit slower, but also
            !gives information about the stability of the HF solution in all directions.

            write(6,"(A)") 
            write(6,"(A)") "Calculating RPA from stability matrix..."

#ifdef __CMPLX
            call stop_all(t_r,"Not coded up for complex integrals. Bug ghb24")
#endif

            !Stability = ( A  B  )
            !            ( B* A* )
            StabilitySize=2*ov_space
            allocate(Stability(StabilitySize,StabilitySize),stat=ierr)
            Stability(:,:)=0.0_dp
            Stability(1:ov_space,1:ov_space)=A_mat(1:ov_space,1:ov_space)
!Assume all integrals real to start with
!            Stability(ov_space+1:StabilitySize,1:ov_space)=conjg(B(1:ov_space,1:ov_space))
            Stability(ov_space+1:StabilitySize,1:ov_space)=B_mat(1:ov_space,1:ov_space)    
            Stability(1:ov_space,ov_space+1:StabilitySize)=B_mat(1:ov_space,1:ov_space)
!Assume all integrals real to start with
!            Stability(ov_space+1:StabilitySize,ov_space+1:StabilitySize)=conjg(A_mat(1:ov_space,1:ov_space))
            Stability(ov_space+1:StabilitySize,ov_space+1:StabilitySize)=A_mat(1:ov_space,1:ov_space)

!            call writematrix(Stability,'Stability matrix')

            !Now diagonalize
            !Find optimal space
            allocate(StabilityCopy(StabilitySize,StabilitySize))
            StabilityCopy(:,:)=Stability(:,:)
            allocate(W(StabilitySize),stat=ierr)    !Eigenvalues of stability matrix
            allocate(Work(1))
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            W(:)=0.0_dp
            lWork=-1
            info=0
            call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lwork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)

            do i=1,StabilitySize
                if(W(i).lt.0.0_dp) then
                    write(6,*) i,W(i)
                    call stop_all(t_r,"HF solution not stable. Not local minimum. Recompute HF.")
                endif
            enddo
            if(.not.tdirectRPA) then
                write(6,"(A)") "Stability matrix positive definite. HF solution is minimum. RPA stable"
            endif

            !Now compute S^(1/2), and transform into original basis
            allocate(S_half(StabilitySize,StabilitySize))
            S_half(:,:)=0.0_dp
            do i=1,StabilitySize
                S_half(i,i)=sqrt(W(i))
            enddo
            allocate(temp(StabilitySize,StabilitySize))
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,    &
                S_half,StabilitySize,0.0_dp,temp,StabilitySize)
            call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp,StabilitySize,Stability,   &
                StabilitySize,0.0_dp,S_half,StabilitySize)
            !S_half is now S^1/2 in the original basis

            !Check this by squaring it.
            call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,S_half,    &
                StabilitySize,0.0_dp,temp,StabilitySize)
            do i=1,StabilitySize
                do j=1,StabilitySize
                    if(abs(StabilityCopy(i,j)-temp(i,j)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'S^1/2 not calculated correctly in original basis')
                    endif
                enddo
            enddo

            temp(:,:)=0.0_dp
            do i=1,ov_space
                temp(i,i)=1.0_dp
            enddo
            do i=ov_space+1,StabilitySize
                temp(i,i)=-1.0_dp
            enddo

            allocate(temp2(StabilitySize,StabilitySize))
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
                StabilitySize,0.0_dp,temp2,StabilitySize)
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,S_half, &
                StabilitySize,0.0_dp,temp,StabilitySize)
            !Now diagonalize temp = S^(1/2) (1 0 \\ 0 -1 ) S^(1/2)

            lWork=-1
            allocate(W2(StabilitySize))
            allocate(Work(1))
            W2(:)=0.0_dp
            call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lwork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)
!            call writevector(W2,'Excitation energies')
            ! temp now holds the eigenvectors X~ Y~
            ! W2 runs over StabilitySize eigenvalues (ov_space*2). Therefor we expect redundant pairs of +-W2, corresponding
            ! to pairs of eigenvectors (X^v Y^v) and (X^v* Y^v*) (Same in real spaces).
            do i=1,ov_space
                !This they are listed in order of increasing eigenvalue, we should be able to easily check that they pair up
                if(abs(W2(i)+W2(StabilitySize-i+1)).gt.1.0e-7_dp) then
                    write(6,*) i,StabilitySize-i+1, W2(i), W2(StabilitySize-i+1), abs(W2(i)-W2(StabilitySize-i+1))
                    call stop_all(t_r,"Excitation energy eigenvalues do not pair")
                endif
            enddo

            !We actually have everything we need for the energy already now. However, calculate X and Y too.
            !Now construct (X Y) = S^(-1/2) (X~ Y~)
            !First get S^(-1/2) in the original basis
            S_half(:,:)=0.0_dp
            do i=1,StabilitySize
                S_half(i,i)=-sqrt(W(i))
            enddo
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,S_half, &
                StabilitySize,0.0_dp,temp2,StabilitySize)
            call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,Stability,  &
                StabilitySize,0.0_dp,S_half,StabilitySize)
            !S_half is now S^(-1/2) in the original basis

            !Now multiply S^(-1/2) (X~ y~)
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
                StabilitySize,0.0_dp,temp2,StabilitySize)

            !Check that eigenvectors are also paired.
            !Rotations among degenerate sets will screw this up though
!            do i=1,ov_space
!                write(6,*) "Eigenvectors: ",i,StabilitySize-i+1,W2(i),W2(StabilitySize-i+1)
!                do j=1,StabilitySize
!                    write(6,*) j,temp2(j,i),temp2(j,StabilitySize-i+1)
!                enddo
!            enddo
            
!            call writematrix(temp2,'X Y // Y X',.true.)
            !temp2 should now be a matrix of (Y X)
!                                            (X Y)
!           This is the other way round to normal, but due to the fact that our eigenvalues are ordered -ve -> +ve
!           TODO: Are the signs of this matrix correct?
            allocate(X_stab(ov_space,ov_space)) !First index is (m,i) compound index. Second is the eigenvector index.
            allocate(Y_stab(ov_space,ov_space))
            X_stab(:,:)=0.0_dp
            Y_stab(:,:)=0.0_dp
            !Put the eigenvectors corresponding to *positive* eigenvalues into the X_stab and Y_stab arrays.
            X_stab(1:ov_space,1:ov_space)=temp2(1:ov_space,ov_space+1:StabilitySize)
            Y_stab(1:ov_space,1:ov_space)=-temp2(ov_space+1:StabilitySize,ov_space+1:StabilitySize)
            deallocate(temp2)

            !Normalize the eigenvectors appropriately
            do mu=1,ov_space
                norm=0.0_dp
                Y_norm = 0.0_dp
                X_norm = 0.0_dp
                do i=1,ov_space
                    norm = norm + X_stab(i,mu)*X_stab(i,mu) - Y_stab(i,mu)*Y_stab(i,mu)
                    Y_norm = Y_norm + Y_stab(i,mu)*Y_stab(i,mu)
                    X_norm = X_norm + X_stab(i,mu)*X_stab(i,mu)
                enddo
                if(norm.le.0.0_dp) then 
                    write(6,*) "Norm^2 for vector ",mu," is: ",norm
                    call stop_all(t_r,'norm undefined')
                endif
                norm = sqrt(norm)
                do i=1,ov_space
                    X_stab(i,mu) = X_stab(i,mu)/norm
                    Y_stab(i,mu) = Y_stab(i,mu)/norm
                enddo
                if(Y_norm.gt.X_norm/2.0_dp) then
                    write(6,*) "Warning: hole amplitudes large for excitation: ",mu,    &
                        " Quasi-boson approximation breaking down."
                    write(6,*) "Norm of X component: ",X_norm
                    write(6,*) "Norm of Y component: ",Y_norm
                endif
            enddo
!            call writematrix(X_stab,'X',.true.)

            !Now check orthogonality 
            call Check_XY_orthogonality(X_stab,Y_stab)
            
!            call writevector(W2,'Stab_eigenvalues')

            !Now check that we satisfy the original RPA equations
            !For the *positive* eigenvalue space (since we have extracted eigenvectors corresponding to this), check that:
!           ( A B ) (X) = E_v(X )
!           ( B A ) (Y)      (-Y)
            deallocate(temp)
            allocate(temp(StabilitySize,ov_space))
            temp(1:ov_space,1:ov_space) = X_stab(:,:)
            temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
            allocate(temp2(StabilitySize,ov_space))
            call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
                StabilitySize,0.0_dp,temp2,StabilitySize)
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                        write(6,*) i,j,temp2(j,i),(W2(i+ov_space)*X_stab(j,i)),W2(i+ov_space)
                        call stop_all(t_r,"RPA equations not satisfied for positive frequencies in X matrix")
                    endif
                enddo
            enddo
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp2(j+ov_space,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                        write(6,*) i,j,temp2(j+ov_space,i),(-W2(i+ov_space)*Y_stab(j,i)),-W2(i+ov_space)
                        call stop_all(t_r,"RPA equations not satisfied for positive frequencies in Y matrix")
                    endif
                enddo
            enddo
            deallocate(temp,temp2)

            !Is is also satisfied the other way around?
            !Check that we also satisfy (still for the *positive* eigenvalues):
!           ( A B ) (Y) = -E_v(Y )
!           ( B A ) (X)       (-X)
            allocate(temp(StabilitySize,ov_space))
            temp(1:ov_space,1:ov_space) = Y_stab(:,:)
            temp(ov_space+1:StabilitySize,1:ov_space) = X_stab(:,:)
            allocate(temp2(StabilitySize,ov_space))
            call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
                StabilitySize,0.0_dp,temp2,StabilitySize)
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp2(j,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                        call stop_all(t_r,"RPA equations not satisfied for negative frequencies in X matrix")
                    endif
                enddo
            enddo
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp2(j+ov_space,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                        call stop_all(t_r,"RPA equations not satisfied for negative frequencies in Y matrix")
                    endif
                enddo
            enddo
            deallocate(temp,temp2)

            !TODO: Finally, check that we satisfy eq. 1 in the Scuseria paper for X and Y defined for positive eigenvalues...
            do i=ov_space+1,StabilitySize
                do j=1,StabilitySize
                    StabilityCopy(i,j)=-StabilityCopy(i,j)
                enddo
            enddo
            !Stability copy is now (A B // -B -A)
            allocate(temp(StabilitySize,ov_space))
            allocate(temp2(StabilitySize,ov_space))
            temp=0.0_dp
            temp(1:ov_space,1:ov_space) = X_stab(:,:)
            temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
            call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
                StabilitySize,0.0_dp,temp2,StabilitySize)
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-7_dp) then
                        call stop_all(t_r,"RPA equations not satisfied for X")
                    endif
                enddo
            enddo
            do i=1,ov_space
                do j=ov_space+1,StabilitySize
                    if(abs(temp2(j,i)-(W2(i+ov_space)*Y_stab(j-ov_space,i))).gt.1.0e-7_dp) then
                        call stop_all(t_r,"RPA equations not satisfied for Y")
                    endif
                enddo
            enddo
            deallocate(temp,temp2)

            !Now calculate energy, in two different ways:
            !1. -1/2 Tr[A] + 1/2 sum_v E_v(positive)
            Energy_stab=0.0_dp
            do i=1,ov_space
                Energy_stab = Energy_stab + W2(ov_space+i) - A_mat(i,i)
            enddo
            Energy_stab = Energy_stab/2.0_dp

            if(tDirectRPA) then
                write(6,"(A,G25.10)") "Direct RPA energy from stability analysis (plasmonic RPA-TDA excitation energies): ", &
                    Energy_stab
            else
                write(6,"(A,G25.10)") "Full RPA energy from stability analysis (plasmonic RPA-TDA excitation energies): ",  &
                    Energy_stab
            endif

            Energy_stab = 0.0_dp
            !E = 0.25 * Tr[BZ] where Z = Y X^-1

            allocate(temp2(ov_space,ov_space))
            temp2(:,:) = 0.0_dp
            !Find X^-1 
            call d_inv(X_stab,temp2)
!            call writematrix(temp2,'X^-1',.true.)
            allocate(temp(ov_space,ov_space))
            !Find Z (temp)
            call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,Y_stab,ov_space,temp2,ov_space,0.0_dp,temp,ov_space)
            !Find BZ (temp2)
            call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,B_mat,ov_space,temp,ov_space,0.0_dp,temp2,ov_space)
            !Take trace of BZ
            do i=1,ov_space
                Energy_stab = Energy_stab + temp2(i,i)
            enddo
            Energy_stab = Energy_stab/2.0_dp
            if(tDirectRPA) then
                write(6,"(A,G25.10)") "Direct RPA energy from stability analysis (Ring-CCD: 1/2 Tr[BZ]): ",Energy_stab
            else
                write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Ring-CCD: 1/2 Tr[BZ]): ",Energy_stab
            endif

            Energy_stab = 0.0_dp
            do i=1,ov_space
                Y_norm = 0.0_dp
                do j=1,ov_space
                    Y_norm = Y_norm + Y_stab(j,i)**2
                enddo
                Energy_stab = Energy_stab - W2(i+ov_space)*Y_norm
            enddo
            if(tDirectRPA) then
                write(6,"(A,G25.10)") "Direct RPA energy from stability analysis (Y-matrix): ",Energy_stab
            else
                write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Y-matrix): ",Energy_stab
            endif
                
            deallocate(W2,W,temp,temp2,StabilityCopy,Stability)

        endif

        write(6,*) 
        write(6,"(A)") "Calculating RPA via Cholesky decomposition"
        !Now, we calculate the RPA amplitudes via consideration of the boson operators.
        !This will generally be quicker, and should give the same result.

        !Calculate A-B
        allocate(AminB(ov_space,ov_space))
        do i=1,ov_space
            do j=1,ov_space
                AminB(j,i) = A_mat(j,i) - B_mat(j,i)
            enddo
        enddo

        if(.true.) then
            !Optional
            !Check that AminB is positive definite - required for Cholesky decomposition
            allocate(temp(ov_space,ov_space))
            temp(:,:) = AminB(:,:)
            !Also check it is symmetric!
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(AminB(i,j)-AminB(j,i)).gt.1.0e-7_dp) then
                        call stop_all(t_r,"A - B matrix is not symmetric")
                    endif
                enddo
            enddo

            !Now diagonalise
            lWork=-1
            allocate(W(ov_space))
            allocate(Work(1))
            W(:)=0.0_dp
            call dsyev('V','U',ov_space,temp,ov_space,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',ov_space,temp,ov_space,W,Work,lWork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)
            do i=1,ov_space
                if(W(i).lt.0.0_dp) then
                    call writevector(W,'AminB eigenvalues')
                    call stop_all(t_r,"A-B not positive definite - cannot do cholesky decomposition")
                endif
            enddo
            deallocate(temp,W)
        endif

        !Construct triangular matrices via cholesky decomposition
        call dpotrf('U',ov_space,AminB,ov_space,info)
        if(info.ne.0) call stop_all(t_r,'Cholesky decomposition failed.')
        !Now set lower triangular part to zero
        do i=1,ov_space
            do j=1,i-1
                AminB(i,j)=0.0_dp
            enddo
        enddo

        if(.true.) then
            !Optional
!            call writematrix(AminB,'T')
            !Check cholesky decomposition successful
            allocate(temp(ov_space,ov_space))
            call dgemm('T','N',ov_space,ov_space,ov_space,1.0_dp,AminB,ov_space,AminB,ov_space,0.0_dp,temp,ov_space)
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(temp(i,j)-(A_mat(i,j)-B_mat(i,j))).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Cholesky decomposition not as expected')
                    endif
                enddo
            enddo
            deallocate(temp)
        endif

        allocate(AplusB(ov_space,ov_space),stat=ierr)
        do i=1,ov_space
            do j=1,ov_space
                AplusB(j,i) = A_mat(j,i) + B_mat(j,i)
            enddo
        enddo

        allocate(temp(ov_space,ov_space))
        call dgemm('N','N',ov_space,ov_space,ov_space,1.0_dp,AminB,ov_space,AplusB,ov_space,0.0_dp,temp,ov_space)
        call dgemm('N','T',ov_space,ov_space,ov_space,1.0_dp,temp,ov_space,AminB,ov_space,0.0_dp,AplusB,ov_space)
        deallocate(temp)
        !Check T (A + B) T^T is actually symmetric
        do i=1,ov_space
            do j=1,ov_space
                if(abs(AplusB(i,j)-AplusB(j,i)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'Not symmetric')
                endif
            enddo
        enddo

        !AplusB is now actually T (A + B) T^T
        !Diagonalize this
        lWork=-1
        allocate(W(ov_space))
        allocate(Work(1))
        W(:)=0.0_dp
        call dsyev('V','U',ov_space,AplusB,ov_space,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',ov_space,AplusB,ov_space,W,Work,lWork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)
        do i=1,ov_space
            if(W(i).lt.0.0_dp) then
                call stop_all(t_r,"Excitation energies^2 not positive. Error in diagonalization.")
            endif
        enddo
!        call writevector(sqrt(W),'Eigenvalues')

        !Excitation energies are now the +- sqrt of eigenvalues
        !Construct W^(-1/2) T^T R   (temp)
        !and W^(1/2) T^-1 R (temp2)
        allocate(temp(ov_space,ov_space))
        call dgemm('T','N',ov_space,ov_space,ov_space,1.0_dp,AminB,ov_space,AplusB,ov_space,0.0_dp,temp,ov_space)
        do i=1,ov_space
            !Mulitply each column by corresponding eigenvalue^(-1/2). Only deal with positive eigenvalues here
            do j=1,ov_space
                temp(j,i) = temp(j,i)*(1.0_dp/(sqrt(sqrt(W(i)))))
            enddo
        enddo

        allocate(temp2(ov_space,ov_space))
        allocate(temp3(ov_space,ov_space))  !Temp3 will hold the inverse of the T matrix temporarily.

        call d_inv(AminB,temp3)
        call dgemm('N','N',ov_space,ov_space,ov_space,1.0_dp,temp3,ov_space,AplusB,ov_space,0.0_dp,temp2,ov_space)
        deallocate(temp3)

        do i=1,ov_space
            !Mulitply each column by corresponding eigenvalue^(1/2). Only deal with positive eigenvalues here
            do j=1,ov_space
                temp2(j,i) = temp2(j,i)*(sqrt(sqrt(W(i))))
            enddo
        enddo

        allocate(X_Chol(ov_space,ov_space))
        allocate(Y_Chol(ov_space,ov_space))
        do i=1,ov_space
            do j=1,ov_space
                X_Chol(j,i) = 0.5_dp * (temp(j,i) + temp2(j,i))
                Y_Chol(j,i) = 0.5_dp * (temp(j,i) - temp2(j,i))
            enddo
        enddo
!        do i=1,ov_space
!            do j=1,ov_space
!                if(abs(Y_Chol(i,j)).gt.1.0e-7) then
!                    write(6,*) "Found non-zero Y matrix value..."
!                endif
!            enddo
!        enddo

        !Y_Chol and X_Chol now cover all eigenvectors.
        !Eigensystem can be represented as:
        ! ( X ) val=sqrt(W)     and     ( Y* ) val=-sqrt(W)
        ! ( Y )                         ( X* )

        call Check_XY_orthogonality(X_Chol,Y_Chol)

        !Are X and Y vectors the same as from the stability matrix?
        !I guess there will be different rotation among degenerate sets, so not a well defined comparison.
!        if(tStabilityAnalysis) then
!            do i=1,ov_space
!                do j=1,ov_space
!                    if(abs(X_Chol(j,i)-X_Stab(j,i)).gt.1.0e-7) then
!                        write(6,*) j,i,X_Chol(j,i),X_Stab(j,i)
!!                        call stop_all(t_r,"Calculation of X matrix not the same as via stability matrix")
!                    endif
!                    if(abs(Y_Chol(j,i)-Y_Stab(j,i)).gt.1.0e-7) then
!                        write(6,*) j,i,Y_Chol(j,i),Y_Stab(j,i)
!!                        call stop_all(t_r,"Calculation of Y matrix not the same as via stability matrix")
!                    endif
!                enddo
!            enddo
!        endif


        !Now check that we satisfy the original RPA equations
        !For the *positive* eigenvalue space (since we have extracted eigenvectors corresponding to this), check that:
!           ( A B ) (X) = E_v(X )
!           ( B A ) (Y)      (-Y)
        deallocate(temp,temp2)
        allocate(Stability(StabilitySize,StabilitySize))
        Stability(:,:)=0.0_dp
        Stability(1:ov_space,1:ov_space)=A_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,1:ov_space)=B_mat(1:ov_space,1:ov_space)    
        Stability(1:ov_space,ov_space+1:StabilitySize)=B_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,ov_space+1:StabilitySize)=A_mat(1:ov_space,1:ov_space)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = X_Chol(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_Chol(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,Stability,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(sqrt(W(i))*X_Chol(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j,i),(sqrt(W(i))*X_Chol(j,i)),sqrt(W(i))
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(-sqrt(W(i))*Y_Chol(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j+ov_space,i),(-sqrt(W(i))*Y_Chol(j,i)),-sqrt(W(i))
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !Is is also satisfied the other way around?
        !Check that we also satisfy (still for the *positive* eigenvalues):
!           ( A B ) (Y) = -E_v(Y )
!           ( B A ) (X)       (-X)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = Y_Chol(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = X_Chol(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,Stability,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(-sqrt(W(i))*Y_Chol(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(sqrt(W(i))*X_Chol(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !TODO: Finally, check that we satisfy eq. 1 in the Scuseria paper for X and Y defined for positive eigenvalues...
        do i=ov_space+1,StabilitySize
            do j=1,StabilitySize
                Stability(i,j)=-Stability(i,j)
            enddo
        enddo
        !Stability copy is now (A B // -B -A)
        allocate(temp(StabilitySize,ov_space))
        allocate(temp2(StabilitySize,ov_space))
        temp=0.0_dp
        temp(1:ov_space,1:ov_space) = X_Chol(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_Chol(:,:)
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,Stability,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(sqrt(W(i))*X_Chol(j,i))).gt.1.0e-7_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for X")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=ov_space+1,StabilitySize
                if(abs(temp2(j,i)-(sqrt(W(i))*Y_Chol(j-ov_space,i))).gt.1.0e-7_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for Y")
                endif
            enddo
        enddo
        deallocate(temp,temp2,Stability)

        Energy2 = real(HDiagTemp,dp)
        norm = 0.0_dp
        do i=1,ov_space
            norm = norm + A_mat(i,i)
        enddo
        Energy2 = Energy2 - norm/2.0_dp

        norm = 0.0_dp
        do i=1,ov_space
            !Run over positive eigenvalue pairs
            norm = norm + sqrt(W(i))
        enddo
        Energy2 = Energy2 + norm/2.0_dp
        if(tDirectRPA) then
            write(6,"(A,G25.10)") "Direct RPA energy (plasmonic RPA-TDA excitation energies): ",    &
                Energy2-real(HDiagTemp,dp)
        else
            write(6,"(A,G25.10)") "Full RPA energy (plasmonic RPA-TDA excitation energies): ",  &
                Energy2-real(HDiagTemp,dp)
        endif
        
        Energy = 0.0_dp
        allocate(temp(ov_space,ov_space))
        allocate(temp2(ov_space,ov_space))
        temp(:,:) = 0.0_dp
        call d_inv(X_Chol,temp)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,Y_Chol,ov_space,temp,ov_space,0.0_dp,temp2,ov_space)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,B_Mat,ov_space,temp2,ov_space,0.0_dp,temp,ov_space)
        do i=1,ov_space
            Energy = Energy + temp(i,i)
        enddo
        Energy = Energy/2.0_dp
        if(tDirectRPA) then
            write(6,"(A,G25.10)") "Direct RPA energy (Ring-CCD: 1/2 Tr[BZ]): ",Energy
        else
            write(6,"(A,G25.10)") "Full RPA energy (Ring-CCD: 1/2 Tr[BZ]): ",Energy
        endif
        deallocate(temp,temp2)
        

        HDiagTemp = get_helement(fDet, fDet, 0)
        Energy = real(HDiagTemp,dp)
        do v=1,ov_space
            norm=0.0_dp
            do i=1,ov_space
                norm = norm + (abs(Y_Chol(i,v)))**2.0_dp
            enddo
            Energy = Energy - norm*sqrt(W(v))
        enddo
        if(tDirectRPA) then
            write(6,"(A,G25.10)") "Direct RPA energy (Y matrix): ",Energy-real(HDiagTemp,dp)
        else
            write(6,"(A,G25.10)") "Full RPA energy (Y matrix): ",Energy-real(HDiagTemp,dp)
        endif

        write(6,"(A)")
        write(6,"(A)") "RPA calculation completed successfully"
        write(6,"(A)")

        write(6,"(A)") "Calculating 1RDM..."
        allocate(OccNumbers(nBasis))    !This is the GS 1RDM, which is diagonal in RPA

        OccNumbers(:) = 0.0_dp

        !For occupied orbital: Gamma(i,i) = 1 - 1/2 \sum_{mu,a} |Y^{mu}_ia|^2
        !For virtual orbital:  Gamma(a,a) =     1/2 \sum_{mu,i} |Y^{mu}_ia|^2

        do i=1,nel
            do mu=1,ov_space
                do a=virt_start,nBasis
                    ia_ind = ov_space_ind(a,i)
                    OccNumbers(i) = OccNumbers(i) + 0.5_dp*abs(Y_Chol(ia_ind,mu))**2.0_dp
                enddo
            enddo
            OccNumbers(i)=1.0_dp - OccNumbers(i)
        enddo

        do a=virt_start,nBasis
            do mu=1,ov_space
                do i=1,nel
                    ia_ind = ov_space_ind(a,i)
                    OccNumbers(a) = OccNumbers(a) + 0.5_dp*abs(Y_Chol(ia_ind,mu))**2.0_dp
                enddo
            enddo
        enddo

        write(6,"(A)") "Writing RPA occupation numbers to file: RPA_OccNumbers"
        write(6,"(A)") ""
        iunit = get_free_unit()
        open(iunit,file='RPA_OccNumbers',status='unknown')
        write(iunit,"(A)") "# Orbital(Energy_order)  Orbital   Fock eigenvalue   RPA_Occupation"
        do i=1,nBasis
            i_ind = Brr(i)
            write(iunit,"(2I9,2G25.10)") i,i_ind,Arr(i_ind,2),OccNumbers(i)
        enddo
        close(iunit)

        deallocate(A_mat,B_mat,X_Chol,Y_Chol,OccNumbers)

    end subroutine RunRPA_QBA

    SUBROUTINE d_inv (mat,matinv)
        implicit none
        real(dp), INTENT(IN) :: mat(:,:)
        real(dp), dimension(size(mat,1),size(mat,2)), intent(out) :: matinv
        real(dp), dimension(size(mat,1),size(mat,2)) :: matdum
        integer, dimension(size(mat,1)) :: ipiv
        integer :: nsize,msize,i,info

        msize=size(mat,1)
        nsize=size(mat,2)
        matdum=mat
        matinv=0.0_dp
        do i=1,msize
            matinv(i,i)=1.0_dp
        enddo  
        info=0 
        call dGETRF(msize,nsize,matdum,nsize,ipiv,info)
!        IF (INFO /= 0) STOP 'Error with d_inv 1'
        if(info /= 0) then
            write(6,*) 'Warning from d_inv 1', info
            call stop_all("d_inv","Warning from d_inv 1")
        endif 
        call dGETRS('n',msize,nsize,matdum,nsize,IPIV,matinv,msize,info)
        if(info.ne.0) call stop_all("d_inv",'Error with d_inv 2')

    END SUBROUTINE d_inv


    subroutine Check_XY_orthogonality(X,Y)
        implicit none
        real(dp), intent(in) :: X(ov_space,ov_space),Y(ov_space,ov_space)
        integer :: i,i_p,m,m_p,v,mi_ind,mp_ip_ind,mu,mu_p,j
        real(dp) :: temp
        character(len=*), parameter :: t_r="Check_XY_orthogonality"

        !Eigenvectors orthonormal
        do mu=1,ov_space
            do mu_p=1,ov_space
                temp = 0.0_dp
                do i=1,ov_space
                    temp = temp + X(i,mu)*X(i,mu_p) - Y(i,mu)*Y(i,mu_p)
                enddo
                if((mu.eq.mu_p).and.((abs(temp)-1.0_dp).gt.1.0e-6_dp)) then
                    write(6,*) mu,mu_p,temp
                    call writevector(X(:,mu),'X(:,mu)')
                    call writevector(Y(:,mu),'Y(:,mu)')
                    call stop_all(t_r,'X/Y not normalized')
                elseif((mu.ne.mu_p).and.(abs(temp).gt.1.0e-6_dp)) then
                    write(6,*) mu,mu_p
                    call stop_all(t_r,'X/Y not orthogonal')
                endif
            enddo
        enddo

        !Rows orthonormal too
        do i=1,ov_space
            do j=1,ov_space
                temp = 0.0_dp
                do mu=1,ov_space
                    temp = temp + X(i,mu)*X(j,mu) - Y(i,mu)*Y(j,mu)
                enddo
                if((i.ne.j).and.(abs(temp).gt.1.0e-7_dp)) then
                    write(6,*) i,j,temp
                    call stop_all(t_r,'X/Y rows not orthogonal')
                elseif((i.eq.j).and.(abs(temp)-1.0_dp).gt.1.0e-7_dp) then
                    write(6,*) i,j,temp
                    call stop_all(t_r,'X/Y rows not normalized')
                endif
            enddo
        enddo

    end subroutine Check_XY_orthogonality

    !a is fast
    integer function ov_space_ind(a,i)
        implicit none
        integer, intent(in) :: i,a

        ov_space_ind = (i-1)*(nBasis-NEl) + a - NEl

    end function ov_space_ind

    subroutine WriteMatrix(mat,matname,tOneLine)
        implicit none
        real(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(G25.10)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,G25.10)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix

    subroutine WriteVector(vec,vecname)
        implicit none
        real(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector


end module RPA_Mod
