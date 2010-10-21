module ClawData

  type griddat
      double precision, allocatable :: amdq(:,:),apdq(:,:),dtdx(:)
      double precision, allocatable :: amdq2(:,:),apdq2(:,:)
      double precision, allocatable :: qr(:,:),ql(:,:)
      double precision, allocatable :: s(:,:),wave(:,:,:)
      double precision, allocatable :: evl(:,:,:),evr(:,:,:)
      double precision, allocatable :: auxl(:,:),auxr(:,:)
      !Need to do this a better way:
      double precision, allocatable :: dq(:,:),u(:,:,:),hh(:,:),uu(:,:)
      double precision, allocatable :: gg(:,:),uh(:,:,:),dq1m(:)
      ! For 2D:
      double precision, allocatable :: q1d(:,:),dq1d(:,:),aux2(:,:)
      double precision, allocatable :: dtdx1d(:), dtdy1d(:)

      ! For multidimensional reconstruction in 2D:
      double precision, allocatable :: hgg(:,:,:), hff(:,:,:)
      double precision, allocatable :: du(:,:,:),uhh(:,:,:,:),ugg(:,:,:,:)
  end type griddat

end module
