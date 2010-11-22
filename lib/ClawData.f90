module ClawData

  type griddat
      double precision, allocatable :: amdq(:,:),apdq(:,:),dtdx(:)
      double precision, allocatable :: amdq2(:,:),apdq2(:,:)
      double precision, allocatable :: qr(:,:),ql(:,:)
      double precision, allocatable :: s(:,:),wave(:,:,:)
      double precision, allocatable :: evl(:,:,:),evr(:,:,:)
      double precision, allocatable :: auxl(:,:),auxr(:,:)
      ! For 2D:
      double precision, allocatable :: q1d(:,:),dq1d(:,:),aux2(:,:)
      double precision, allocatable :: dtdx1d(:), dtdy1d(:)

      double precision, allocatable :: qgauss(:,:,:,:), q1dgauss(:,:,:)
      double precision, allocatable :: qlgauss(:,:,:), qrgauss(:,:,:)
  end type griddat

end module
