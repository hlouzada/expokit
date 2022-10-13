def dgesl():
    """dgesl solves the double preision system
     a * x = b  or  trans(a) * x = b
     using the fators omputed by dgeo or dgefa.

     on entry

        a       double preision(lda, n)
                the output from dgeo or dgefa.

        lda     integer
                the leading dimension of the array  a .

        n       integer
                the order of the matrix  a .

        ipvt    integer(n)
                the pivot vetor from dgeo or dgefa.

        b       double preision(n)
                the right hand side vetor.

        job     integer
                = 0         to solve  a*x = b ,
                = nonzero   to solve  trans(a)*x = b  where
                            trans(a)  is the transpose.

     on return

        b       the solution vetor  x .

     error ondition

        a division by zero will our if the input fator ontains a
        zero on the diagonal.  tehnially this indiates singularity
        but it is often aused by improper arguments or improper
        setting of lda .  it will not our if the subroutines are
        alled orretly and if dgeo has set rond .gt. 0.0
        or dgefa has set info .eq. 0 .

     to ompute  inverse(a) *   where    is a matrix
     with  p  olumns
           all dgeo(a,lda,n,ipvt,rond,z)
           if (rond is too small) go to ...
           do 10 j = 1, p
              all dgesl(a,lda,n,ipvt,(1,j),0)
        10 ontinue

     linpak. this version dated 08/14/78 .
     leve moler, university of new mexio, argonne national lab.

     subroutines and funtions

     blas daxpy,ddot

     internal variables"""
    ...