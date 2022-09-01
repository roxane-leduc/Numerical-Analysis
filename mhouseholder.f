      program mhouseholder
              double precision, dimension(:,:), allocatable :: A,T
              integer i,j,ichoix

              read(*,*) n

              allocate(A(n,n))
              allocate(T(n,n))

              read(*,*) ichoix
              if (ichoix ==0) then
                do i=1,n
                        read(*,*) A(i,:)
                enddo
              else
                      do i=1,n
                        do j=1,n
                                A(i,j) = 1.d0 / (i + j - 1.d0)
                         enddo
                       enddo
              endif

              write(*,*) "Phase 1 : Triangularisation"
              call triangularization(A,n,T)

              write(*,*) "Matrice T:"
              write(*,*) ((T(i,j), " ", j=1,n), new_line(""), i=1,n )

              write(*,*) "Phase 2 : Valeurs propres de A"
              call eigenvalues(T,n)

              deallocate(A)
              deallocate(T)
      end
