!TODO: Make sure everything has an ierr argument when allocating, and that everything is deallocated at end
!TODO: Fix RPA from stability matrix
!TODO: Add MP2 for comparison
!TODO: Add symmetry information by using singles excitation generators (store excitations)? (Not urgent) 
!Module to non-iteratively calculate the RPA energy under the quasi-boson approximation
module RPA_Mod
    use SystemData, only: nel, nBasis, Arr, Brr
    use sltcnd_mod, only: sltcnd_2
    use constants, only: dp, int64, n_int
    use Determinants, only: get_helement, fDet
    use SymExcit3, only: GenExcitations3
    use Determinants, only: GetH0Element3
    use bit_reps, only: NIfTot
    use DetBitops, only: EncodeBitDet

    implicit none

    !input options
    logical :: tStabilityAnalysis=.false.

    !global variables
    integer :: ov_space
    integer :: virt_start

    !global arrays
    real(dp), allocatable :: A_mat(:,:)
    real(dp), allocatable :: B_mat(:,:)
    
    contains

    subroutine RunRPA_QBA(Weight,Energy)
        implicit none
        integer :: ierr,i,j,m,n,ex(2,2),ex2(2,2),mi_ind,nj_ind
        integer :: StabilitySize,lWork,info,i_p,m_p,v,mp_ip_ind,ic
        integer :: nJ(NEl),exflag
        integer(n_int) :: iLutHF(0:NIfTot)
        real(dp), intent(out) :: Weight,Energy
        real(dp) :: Energy_stab,Temp_real,norm,Energy2,H0tmp,Fii
        HElement_t :: HDiagTemp,hel
        logical :: tAllExcitsFound,tParity
        real(dp), allocatable :: Stability(:,:),temp2(:,:),W2(:)
        real(dp), allocatable :: W(:),Work(:),S_minhalf(:,:),temp(:,:)
        real(dp), allocatable :: X_stab(:,:),Y_stab(:,:)
        real(dp), allocatable :: X_chol(:,:),Y_chol(:,:)
        real(dp), allocatable :: AminB(:,:),AplusB(:,:),temp3(:,:)
        character(len=*), parameter :: t_r="RunRPA_QBA"

        write(6,"(A)")
        write(6,"(A)") "*****************************************************************"
        write(6,"(A)") "*   Entering RPA calculation with quasi-boson approximation...  *"
        write(6,"(A)") "*****************************************************************"
        write(6,"(A)") 
        call neci_flush(6)

        !Quickly find correlation energy from MP2 for comparison.
        Temp_real = 0.0_dp
        tAllExcitsFound=.false.
        exflag=3
        ex(:,:)=0
        call EncodeBitDet(FDet,iLutHF)
        HDiagTemp=GetH0Element3(FDet)
        Fii=real(HDiagTemp,dp)
        do while(.true.)
            call GenExcitations3(FDet,iLutHF,nJ,exflag,Ex,tParity,tAllExcitsFound)
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

                        A_mat(mi_ind,nj_ind) = sltcnd_2(ex,.false.)
                        B_mat(mi_ind,nj_ind) = sltcnd_2(ex2,.false.)

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
#ifdef __CMPLX
!                if(abs(A_mat(i,j)-conjg(A_mat(j,i))).gt.1.0e-7_dp) then
#else
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7_dp) then
#endif
                    write(6,*) i,j,A_mat(i,j), A_mat(j,i),abs(A_mat(i,j)-A_mat(j,i))
                    
                    call stop_all(t_r,"A not hermitian")
                endif
            enddo
        enddo

        if(tStabilityAnalysis) then
            !Calculate stability analysis of HF solution, and from this, calculate RPA amplitudes
            !This should give identical results to other method, will be quite a bit slower, but also
            !gives information about the stability of the HF solution in all directions.

            write(6,"(A)") "Calculating RPA from stability matrix"

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
            allocate(W(StabilitySize),stat=ierr)    !Eigenvalues of stability matrix
            allocate(Work(1))
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            W(:)=0.0_dp
            lWork=-1
            call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=work(1)
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
            write(6,"(A)") "Stability matrix positive definite. HF solution is minimum. RPA stable"

            !Now compute S^(1/2), and transform into original basis
            allocate(S_minhalf(StabilitySize,StabilitySize))
            S_minhalf(:,:)=0.0_dp
            do i=1,StabilitySize
                S_minhalf(i,i)=sqrt(W(i))
            enddo
            allocate(temp(StabilitySize,StabilitySize))
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,S_minhalf,StabilitySize,0.0_dp,temp,StabilitySize)
            call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp,StabilitySize,Stability,StabilitySize,0.0_dp,S_minhalf,StabilitySize)
            !S_minhalf is now S^1/2 in the original basis

            temp(:,:)=0.0_dp
            do i=1,ov_space
                temp(i,i)=1.0_dp
            enddo
            do i=ov_space+1,StabilitySize
                temp(i,i)=-1.0_dp
            enddo

            allocate(temp2(StabilitySize,StabilitySize))
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_minhalf,StabilitySize,temp,StabilitySize,0.0_dp,temp2,StabilitySize)
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,S_minhalf,StabilitySize,0.0_dp,temp,StabilitySize)
            !Now diagonalize temp = S^(1/2) (1 0 \\ 0 -1 ) S^(1/2)

            lWork=-1
            allocate(W2(StabilitySize))
            allocate(Work(1))
            W2(:)=0.0_dp
            call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=work(1)
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
            S_minhalf(:,:)=0.0_dp
            do i=1,StabilitySize
                S_minhalf(i,i)=-sqrt(W(i))
            enddo
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,S_minhalf,StabilitySize,0.0_dp,temp2,StabilitySize)
            call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,Stability,StabilitySize,0.0_dp,S_minhalf,StabilitySize)
            !S_minhalf is now S^(-1/2) in the original basis

            !Now multiply S^(-1/2) (X~ y~)
            call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_minhalf,StabilitySize,temp,StabilitySize,0.0_dp,temp2,StabilitySize)

            !temp2 now runs over X, Y
            allocate(X_stab(ov_space,ov_space)) !First index is (m,i) compound index. Second is the eigenvector index.
            allocate(Y_stab(ov_space,ov_space))
            X_stab(:,:)=0.0_dp
            Y_stab(:,:)=0.0_dp
            !Put these into the X_stab and Y_stab arrays.
            X_stab(1:ov_space,1:ov_space)=temp2(1:ov_space,1:ov_space)
            Y_stab(1:ov_space,1:ov_space)=temp2(ov_space+1:StabilitySize,1:ov_space)
            deallocate(temp2)

!            call writevector(W2,'Stab_eigenvalues')

            !Now check orthogonality (normalization currently commented out)
!            call Check_XY_orthogonality(X_stab,Y_stab)

            !Now calculate energy, in two different ways
!            call CalcRPAEnergy(Y_stab,W2,Energy_stab)

            Energy = Energy_stab

            deallocate(W2,W,temp)

        endif

        write(6,"(A)") "Calculating RPA from consideration of boson operators"
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
            lwork=work(1)
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
                AplusB(j,i) = A_mat(j,i) - B_mat(j,i)
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
        lwork=work(1)
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',ov_space,AplusB,ov_space,W,Work,lWork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)
        do i=1,ov_space
            if(W(i).lt.0.0_dp) then
                call stop_all(t_r,"Excitation energies not positive. Error in diagonalization.")
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

        if(.true.) then
        !Check orthogonality
            do i=1,nel
                do i_p=1,nel
                    do m=virt_start,nBasis
                        do m_p=virt_start,nBasis
                            temp_real=0.0_dp
                            mi_ind = ov_space_ind(m,i)
                            mp_ip_ind = ov_space_ind(m_p,i_p)
                            do v=1,ov_space
                                temp_real=temp_real+(X_Chol(mi_ind,v)*X_Chol(mp_ip_ind,v)-Y_Chol(mi_ind,v)*Y_Chol(mp_ip_ind,v))
                            enddo
                            if((i.eq.i_p).and.(m.eq.m_p)) then
                                !Diagonal component. Should = 1?
                                if(abs(temp_real-1.0_dp).gt.1.0e-7_dp) then
!                                    write(6,*) i,i_p,m,m_p,temp
                                    call stop_all(t_r,"X/Y eigenvectors not normalized")
                                endif
                            else
                                !Orthogonality. Should = 0
                                if(abs(temp_real).gt.1.0e-7_dp) then
                                    write(6,*) i,i_p,m,m_p,temp
                                    call stop_all(t_r,"X/Y eigenvectors not orthogonal")
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        
        HDiagTemp = get_helement(fDet, fDet, 0)
        Energy = real(HDiagTemp,dp)

        do v=1,ov_space
            norm=0.0_dp
            do i=1,ov_space
                norm = norm + (abs(Y_Chol(i,v)))**2.0_dp
            enddo
            Energy = Energy - norm*sqrt(W(v))
        enddo
        write(6,"(A,G25.10)") "RPA correlation energy from Y matrix: ",Energy-real(HDiagTemp,dp)
        
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
        write(6,"(A,G25.10)") "RPA correlation energy from A matrix boson formulation: ",Energy2-real(HDiagTemp,dp)

        write(6,"(A)") "RPA calculation completed successfully"

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
        info=-1
        call dGETRF(msize,nsize,matdum,nsize,ipiv,info)
!        IF (INFO /= 0) STOP 'Error with d_inv 1'
        if(info /= 0) then
            write(6,*) 'Warning from d_inv 1', info
            call stop_all("d_inv","Warning from d_inv 1")
        endif 
        info=-1
        call dGETRS('n',msize,nsize,matdum,nsize,IPIV,matinv,msize,info)
        if(info.ne.0) call stop_all("d_inv",'Error with d_inv 2')

    END SUBROUTINE d_inv


    subroutine CalcRPAEnergy(Y,W,nrg)
        implicit none
        real(dp) , intent(out) :: nrg
        real(dp) , intent(in) :: Y(ov_space,ov_space*2),W(ov_space*2)
        HElement_t :: HDiagTemp
        real(dp) :: temp,nrg2,norm
        integer :: i,v,m,mi_ind
        character(len=*), parameter :: t_r="CalcRPAEnergy"

        !Calculate energy first way...
        !First sum in HF energy
        HDiagTemp = get_helement(fDet, fDet, 0)
        nrg = real(HDiagTemp,dp)

        !Now include trace of A
        temp = 0.0_dp
        do i=1,ov_space
            temp = temp - A_mat(i,i)
        enddo
        nrg = nrg + temp/2.0_dp

        temp = 0.0_dp
        do i=ov_space+1,ov_space*2
            !Run over positive eigenvalue pairs
            temp = temp + W(i)
        enddo
        nrg = nrg + temp/2.0_dp

        !Now calculate the second way and check the same
        nrg2 = real(HDiagTemp,dp)

        do v=1,ov_space*2
            norm=0.0_dp
            do i=1,nel
                do m=virt_start,nBasis
                    mi_ind = ov_space_ind(m,i)
                    norm = norm + (abs(Y(mi_ind,v)))**2.0_dp
                enddo
            enddo
            nrg2 = nrg2 - norm*W(v)
        enddo

        if(abs(nrg-nrg2).gt.1.0e-7_dp) then
            write(6,*) "nrg: ",nrg
            write(6,*) "nrg2: ",nrg2
            call warning_neci(t_r,"Two RPA energies don't match")
        endif

    end subroutine CalcRPAEnergy


    subroutine Check_XY_orthogonality(X,Y)
        implicit none
        real(dp), intent(in) :: X(ov_space,ov_space),Y(ov_space,ov_space)
        integer :: i,i_p,m,m_p,v,mi_ind,mp_ip_ind
        real(dp) :: temp
        character(len=*), parameter :: t_r="Check_XY_orthogonality"

        do i=1,nel
            do i_p=1,nel
                do m=virt_start,nBasis
                    do m_p=virt_start,nBasis
                        temp=0.0_dp
                        mi_ind = ov_space_ind(m,i)
                        mp_ip_ind = ov_space_ind(m_p,i_p)
                        do v=1,ov_space
                            temp=temp+(X(mi_ind,v)*X(mp_ip_ind,v)-Y(mi_ind,v)*Y(mp_ip_ind,v))
                        enddo
                        if((i.eq.i_p).and.(m.eq.m_p)) then
                            !Diagonal component. Should = 1?
                            if(abs(temp-1.0_dp).gt.1.0e-7_dp) then
!                                write(6,*) i,i_p,m,m_p,temp
!                                call stop_all(t_r,"X/Y eigenvectors not normalized")
                            endif
                        else
                            !Orthogonality. Should = 0
                            if(abs(temp).gt.1.0e-7_dp) then
                                write(6,*) i,i_p,m,m_p,temp
                                call stop_all(t_r,"X/Y eigenvectors not orthogonal")
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine Check_XY_orthogonality

    !m is fast
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
