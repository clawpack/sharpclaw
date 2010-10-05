Module Global

    double precision, allocatable :: q(:,:),aux(:,:)
    double precision, allocatable :: qold(:,:)

  type rkreg
      double precision, allocatable :: qrk(:,:)
  end type rkreg

end module
