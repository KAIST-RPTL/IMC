module rgbimage_m
 
  implicit none
 
  
  public :: rgbimage
 
  type rgbimage
    !! usage
    !!    1) init
    !!    2a) fill_image
    !!      or
    !!    2b) set_pixel
    !!    3) normalize
    !!    4) write
 
    private
 
    integer, dimension(:,:,:), allocatable :: rgb
      !! pixel arrays of rgb values
      !! indices (i,j,k)
      !!    i: position x_i
      !!    j: position y_j
      !!    k=1: red, k=2: green, k=3: blue
 
    integer :: n(2) = 0
      !! image dimensions: [height, width]
	
 
 
  contains
    procedure :: init       => rgbimage_init              ! inits image
    procedure :: fill_image => rgbimage_fill_image        ! fill image with constant rgb value
    procedure :: get_pixel  => rgbimage_get_pixel         ! gets one pixel
    procedure :: normalize  => rgbimage_normalize         ! normalizes all pixels onto range [0, 255]
    procedure :: set_pixel  => rgbimage_set_pixel         ! sets one pixel
    procedure :: write      => rgbimage_write             ! outputs image to file
 
    procedure, private :: inside => rgbimage_inside
    procedure, private :: valid  => rgbimage_valid
  end type
 
 
 	logical :: plotgeom = .false. 
	real(8), allocatable, dimension(:) :: plt_x0, plt_y0, plt_z0, plt_x1, plt_y1, plt_z1
	real(8), allocatable, dimension(:) :: plt_dx, plt_dy, plt_dz
	integer, allocatable, dimension(:) :: plt_nx, plt_ny, plt_nz
	integer :: n_plot = 0
	character(40), allocatable :: plotlist(:)
	integer, allocatable :: plottype(:)
 
 
contains
 
  subroutine rgbimage_init(this, height, width)
    !! initialize image.
    !! sets dimensions, allocates pixels and sets colors to 0.
 
    class(rgbimage), intent(out) :: this
    integer,         intent(in)  :: height, width
 
    this%n = [height, width]
    allocate (this%rgb(3, height,width), source=0)
  end subroutine
 
  logical function rgbimage_valid(this, check_rgb_vals)
    !! checks if the image has valid dimensions and optionally valid rgb values.
 
    class(rgbimage), intent(in)           :: this
    logical,         intent(in), optional :: check_rgb_vals
      !! check if rgb values are in allowed range [0, 255]?
      !! default: dont check
 
    ! always check that dimensions match
    rgbimage_valid = ( all(this%n > 0)                     .and. &
      &               (size(this%rgb, dim=2) == this%n(1)) .and. &
      &               (size(this%rgb, dim=3) == this%n(2)) .and. &
      &               (size(this%rgb, dim=1) == 3)               )
 
    ! optionally: check if rgb values are in allowed range
    if (present(check_rgb_vals)) then
      if (check_rgb_vals) rgbimage_valid = ( rgbimage_valid       .and. &
        &                                   (all(this%rgb >= 0))  .and. &
        &                                   (all(this%rgb <= 255))      )
    end if
 
 
  end function
 
  logical function rgbimage_inside(this, x, y)
    !! checks if given coordinates are inside the image
 
    class(rgbimage), intent(in) :: this
    integer,         intent(in) :: x, y
 
    rgbimage_inside = ((x > 0) .and. (x <= this%n(1)) .and. (y > 0) .and. (y <= this%n(2)))
  end function
 
  subroutine rgbimage_set_pixel(this, x, y, rgb)
    class(rgbimage), intent(inout) :: this
    integer,         intent(in)    :: x, y
      !! coordinates
    integer,         intent(in)    :: rgb(3)
      !! red, green, blue values
 
    if (this%inside(x, y)) then
      ! use given data at first
      this%rgb(:,x,y) = rgb
 
      ! check if given data was out of bounds
      where     (this%rgb(:,x,y) > 255)
        this%rgb(:,x,y) = 255
      elsewhere (this%rgb(:,x,y) < 0)
        this%rgb(:,x,y) = 0
      end where
    end if
  end subroutine
 
  function rgbimage_get_pixel(this, x, y) result(rgb)
    class(rgbimage), intent(in) :: this
    integer,         intent(in) :: x, y
      !! coordinates
    integer                     :: rgb(3)
      !! red, green, blue values
 
    if (this%inside(x, y)) then
      rgb = this%rgb(:,x,y)
    else
      rgb = 0
    end if
  end function
 
  subroutine rgbimage_normalize(this)
    !! normalize colors to be in range [0, 255]
 
    class(rgbimage), intent(inout) :: this
 
    where     (this%rgb(:,:,:) > 255)
      this%rgb(:,:,:) = 255
    elsewhere (this%rgb(:,:,:) < 0)
      this%rgb(:,:,:) = 0
    end where
  end subroutine
 
  subroutine rgbimage_fill_image(this, rgb)
    !! fill whole image with given rgb values.
 
    class(rgbimage), intent(inout) :: this
    integer,         intent(in)    :: rgb(3)
      !! red, green, blue values
 
    integer :: i, j, k
 
    if (this%valid()) then
      do i = 1, 3
        do j = 1, size(this % rgb, dim=2)
            do k = 1, size(this % rgb, dim=3)
                this%rgb(i,j,k) = rgb(i)
            enddo
        enddo
      end do
    end if
  end subroutine
 
  subroutine rgbimage_write(this, fname)
    class(rgbimage), intent(in) :: this
    character(*),    intent(in) :: fname
      !! file path, e.g. "tmp/out.ppm"
 
    integer :: iounit, ios, i,j,k
 
    open (newunit=iounit, file=fname, iostat=ios, action='WRITE')
    if (ios /= 0) error stop "Error opening file: " !// fname
 
    ! write header
    write (iounit, '(A2)')         'P6'
    write (iounit, '(I0, A, I0)') this%n(1), " ", this%n(2)
    write (iounit, '(I0)')         255
 
	do j = this%n(2), 1, -1
		do i = 1, this%n(1)
			write (iounit, '(3A1)', advance='no') this%rgb(:,i,j) 
		end do
	end do
 
 
	!> Original but wrong
    !do i = 1, this%n(1)
    !  do j = 1, this%n(2)
    !    write (iounit, '(3A1)', advance='no') [(achar(this%rgb(i,j,k)), k=1,3)]
    !  end do
    !end do
 
    close (unit=iounit, iostat=ios)
    if (ios /= 0) error stop "Error closing file"
  end subroutine
  
  
  
  
  

  
  
end module
