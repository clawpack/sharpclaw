Module Global

    double precision, allocatable :: q(:,:,:)
    double precision, allocatable :: qold(:,:,:)
    double precision, allocatable :: aux(:,:,:)

    type rkreg
        double precision, allocatable :: qrk(:,:,:)
    end type rkreg

end module
