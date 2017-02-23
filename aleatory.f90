        MODULE aleatory
        REAL, PARAMETER :: pi = 4. *ATAN(1. )
        CONTAINS

        REAL FUNCTION normal(mean,sigma)

          IMPLICIT NONE
          REAL :: mean, sigma, u1, u2

          CALL RANDOM_NUMBER(u1)
          CALL RANDOM_NUMBER(u2)

          normal = sigma*SQRT(-2. *LOG(u1))*COS(2. *pi*u2) + mean

          RETURN
        END FUNCTION normal

        FUNCTION mat_normal(m,n,mean,sigma) RESULT(A)
          IMPLICIT NONE
          INTEGER :: m,n
          REAL :: mean,sigma,A(m,n)
          !local
          REAL :: u1,u2
          integer :: i,j
          !comeca rotina

          do i = 1,m
            do j = 1,n
11            CALL RANDOM_NUMBER(u1)
              CALL RANDOM_NUMBER(u2)
              A(i,j) = sigma*SQRT(-2. *LOG(u1))*COS(2. *pi*u2) + mean
              if ((A(i,j) .gt. mean + 3*sigma) .or. (A(i,j) .lt. mean - 3*sigma)) goto 11
            end do
          end do
        END FUNCTION

        SUBROUTINE seed_from_time(seed)

          IMPLICIT NONE
          INTEGER :: seed, value(8)
          CHARACTER(len=10) :: b(3)
          CALL DATE_AND_TIME( b(1), b(2), b(3), value )
          seed = value(1)+value(2)+value(3)+value(4)+value(5)+value(6)+value(7)+value(8)
          seed = seed + value(1)+value(2)+value(3)+value(4)+value(5)/100+value(6)*100+value(7)/10+value(8)*10

        END SUBROUTINE seed_from_time

        SUBROUTINE init_random_number(iseed)
          INTEGER :: i, seed(12), iseed
          DO i = 1, 12
            seed(i) = i*iseed
          END DO
          CALL RANDOM_SEED(put=seed)
          RETURN
        END SUBROUTINE init_random_number

        INTEGER FUNCTION  rand_int(i,f)
          IMPLICIT NONE
          INTEGER :: i,f
          REAL :: u
          CALL RANDOM_NUMBER(u)
          rand_int = i + FLOOR((f+1-i)*u)  ! We want to choose one from m-n+1 integers
        END FUNCTION

        END MODULE
