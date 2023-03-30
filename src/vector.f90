module vector
  use global_var, only: wp, stdout, PI
  implicit none

  private
  public :: VecR, Mol, VecI, Allocate_Vec, sum, SumVec2

  type VecR
    real(wp), public :: x  !< 向量的x分量
    real(wp), public :: y  !< 向量的y分量
    real(wp), public :: z  !< 向量的y分量
  contains
    !! 任何一个派生类与任意一个不属于该类的数据类型做四则运算, 都需要设置至少2个函数, 满足相关运算法则
    procedure, private, pass(self) :: vec_add_vec, vec_add_real, real_add_vec
    procedure, private, pass(self) :: vec_sub_vec, vec_sub_real, real_sub_vec
    procedure, private, pass(self) :: vec_pmul_vec, vec_mul_real, real_mult_vec
    procedure, private, pass(self) :: vec_equal_vec, vec_equal_real
    procedure, private, pass(self) ::vecr_Div_vecI, vecr_Div_vecr
    procedure, public, pass(self) :: VSAdd, VLen, fprintf, VMul, VRand, VProd, VLensq

    generic :: operator(+) => vec_add_vec, vec_add_real, real_add_vec
    generic :: operator(-) => vec_sub_vec, vec_sub_real, real_sub_vec
    generic :: operator(*) => vec_pmul_vec, vec_mul_real, real_mult_vec
    generic :: assignment(=) => vec_equal_vec, vec_equal_real
    generic :: VDiv => vecr_Div_vecI, vecr_Div_vecr   !! 写到这里的时候才发现变量类型错误, 然后就紧急写了个重载弥补
  end type VecR

  interface sum
    procedure ::SumVectors, SumVector
  end interface

  type Mol
    type(VecR), public :: r      !< 位置向量
    type(VecR), public  :: rv    !< 速度向量
    type(VecR), public  :: ra    !< 加速度向量
  end type Mol

  type VecI
    integer, public :: x
    integer, public :: y
    integer, public :: z
  contains
    procedure, public, pass(self) :: VIProd
  end type

  type Prop
    real(wp) :: val       !< 物理量的值
    real(wp) :: sum       !< 一个物理量的总和, 用于评估平均值
    real(wp) :: sum2      !< 物理量的平方和, 用于评估标准差
  contains
    procedure, public, pass(self) :: PropZero
    procedure, public, pass(self) :: PropAccum
    procedure, public, pass(self) :: PropAvg
  end type

  !! 逐元函数(elemental), 意味着传入的即使是数组也可以进行相关的运算操作
contains

  subroutine Allocate_Vec(Atom, nMol)
    !! 为派生类型数组分配对应的内存空间
    type(Mol), intent(inout), allocatable :: Atom(:)      !< 传入的原子集合信息, 包括位置, 速度, 加速度
    integer, intent(in) ::  nMol                          !< 原子个数
    integer ::  Alloc_flag
    character(len=256) :: Err_message

    if (nMol <= 0) then
      write(stdout, '(A)') "在分配内存时出现了负数的元素个数!"
      stop
    end if

    if (allocated(Atom)) THEN
      deallocate(Atom, stat=Alloc_flag, ERRMSG=Err_message)
      if (Alloc_flag /= 0) then
        write(stdout, '(A)') '在释放Atom数组内存时出错, 错误信息如下: '
        write(stdout, '(A)') Err_message
        stop
      end if
    END IF

    allocate(Atom(nMol), stat=Alloc_flag, ERRMSG=Err_message)
    if (Alloc_flag /= 0) then
      write(stdout, '(A)') '在分配Atom数组内存时出错, 错误信息如下: '
      write(stdout, '(A)') Err_message
      stop
    end if

  end subroutine Allocate_Vec

  pure elemental function vec_add_vec(self, vec2) result(res)
    !! 矢量加法, 两个矢量相加
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    type(VecR) :: res !< 返回结果同样为矢量

    res%x = self%x + vec2%x
    res%y = self%y + vec2%y
    res%z = self%z + vec2%z

  end function vec_add_vec

  pure elemental function vec_add_real(self, var) result(res)
    !! 矢量加法, 将一个矢量与一个实数相加
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x + var
    res%y = self%y + var
    res%z = self%z + var

  end function vec_add_real

  pure elemental function real_add_vec(var,self) result(res)
    !! 矢量加法, 将一个矢量与���个实数相加
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x + var
    res%y = self%y + var
    res%z = self%z + var

  end function real_add_vec

  pure elemental function vec_sub_vec(self, vec2) result(res)
    !! 矢量减法, 两个矢量相减
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    type(VecR) :: res !< 返回结果同样为矢量

    res%x = self%x - vec2%x
    res%y = self%y - vec2%y
    res%z = self%z - vec2%z

  end function vec_sub_vec

  pure elemental function vec_sub_real(self, var) result(res)
    !! 矢量减法, 将一个矢量与一个实数相减
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x - var
    res%y = self%y - var
    res%z = self%z - var

  end  function vec_sub_real

  pure elemental function real_sub_vec(var, self) result(res)
    !! 矢量减法, 将一个矢量与一个实数相减
    class(VecR), intent(in) :: self  !< 输入1个矢量类
    real(wp), intent(in) :: var     !< 要减去的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x - var
    res%y = self%y - var
    res%z = self%z - var

  end function real_sub_vec

  pure elemental function vec_pmul_vec(self, vec2) result(res)
    !! 矢量的点乘运算法则
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    real(wp) :: res !< 返回双精度的结果

    res = self%x * vec2%x + self%y * vec2%y + self%z * vec2%z

  end function vec_pmul_vec

  pure elemental function vec_mul_real(self, var) result(res)
    !! 矢量乘以一个实数
    class(VecR), intent(in) :: self  !< 输入的2个矢量类
    real(wp), intent(in) :: var     !< 要乘以的一个实数
    type(VecR) :: res !< 返回1个缩放var倍的向量

    res%x = self%x * var
    res%y = self%y * var
    res%z = self%z * var

  end function vec_mul_real

  pure elemental function real_mult_vec(var, self) result(res)
    !! 矢量乘以一个实数
    class(VecR), intent(in) :: self  !< 输入的2个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res !< 返回一个缩放之后的矢量

    res%x = self%x * var
    res%y = self%y * var
    res%z = self%z * var

  end function real_mult_vec

  pure elemental subroutine vec_equal_vec(self, vec)
    !! 矢量等于一个矢量
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    type(VecR), intent(in) :: vec         !< 另一个矢量

    self%x = vec%x
    self%y = vec%y
    self%z = vec%z

  end subroutine vec_equal_vec

  pure elemental subroutine vec_equal_real(self, var)
    !! 让矢量等于一个实数
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    real(wp), intent(in) :: var           !< 实数

    self%x = var
    self%y = var
    self%z = var

  end subroutine vec_equal_real

  pure subroutine VSAdd(self, vec, var)
    !! 将self向量加上一个放大了var倍的vec向量
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    type(VecR), intent(in) :: vec         !< 被放大的矢量
    real(wp), intent(in) :: var           !< 放大倍数

    self%x = self%x + vec%x * var
    self%y = self%y + vec%y * var
    self%z = self%z + vec%z * var

  end subroutine VSAdd

  pure elemental function VLen(self) result(res)
    !! 返回矢量的长度
    class(VecR), intent(in) :: self   !< 输入矢量
    real(wp) :: res                   !< 矢量长度

    res = sqrt(self%x * self%x + self%y * self%y + self%z * self%z)

  end function VLen

  pure elemental function VLensq(self) result(res)
    !! 返回矢量的长度的平方
    class(VecR), intent(in) :: self   !< 输入矢量
    real(wp) :: res                   !< 矢量长度

    res = self%x * self%x + self%y * self%y + self%z * self%z

  end function VLensq

  pure elemental subroutine VMul(self, vec1, vec2)
    !! 将vec1�����vec2对应的x和y相乘并赋值给self
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1, vec2

    self%x = vec1%x * vec2%x
    self%y = vec1%y * vec2%y
    self%z = vec1%z * vec2%z

  end  subroutine VMul

  pure elemental subroutine vecr_Div_vecr(self, vec1, vec2)
    !! 将vec1和vec2对应的x和y相除并���值给self
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1, vec2

    self%x = vec1%x / vec2%x
    self%y = vec1%y / vec2%y
    self%z = vec1%z / vec2%z

  end  subroutine vecr_Div_vecr

  pure elemental subroutine vecr_Div_vecI(self, vec1, vec2)
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1
    type(vecI), intent(in) :: vec2

    self%x = vec1%x / vec2%x
    self%y = vec1%y / vec2%y
    self%z = vec1%z / vec2%z

  end subroutine vecr_Div_vecI

  impure elemental subroutine VRand(self)
    !! 将vecr类变量中的分量随机赋予一个值
    class(VecR), intent(inout) :: self
    real(wp) :: seed1, seed2, s
    real(wp) :: x, y, z

    s = 2.d0
    do
      call random_seed()
      call random_number(seed1)
      call random_number(seed2)
      x = 2 * seed1 - 1.d0
      y = 2 * seed2 - 1.d0
      s = x * x + y * y
      if (s < 1.d0) exit
    end do
    self%z =  1.d0 - 2.d0 * s
    seed2 = 2.d0 * sqrt(1.d0 - s)
    self%x =  s * x
    self%y =  s * y

  end subroutine VRand

  pure elemental function VProd(self) result(res)
    !! 返回向量坐标相乘的结果
    class(VecR), intent(in) :: self
    real(wp) :: res

    res = self%x * self%y * self%y
  end function VProd

  subroutine fprintf(self, fileid)
    !! 将向量以格式化的形式输出
    class(VecR), intent(in) :: self
    integer, intent(in), optional :: fileid     !< 可选的通道id, 如果传入的话就可以输出到文件中

    if (.not. present(fileid) ) then
      write(stdout, '(A, F0.7, A, F0.7, A, F0.7, A)') '( ', self%x, ', ', self%y, ', ', self%z, ' )'
    else
      write(fileid, '( I0, 5X, F0.7, 5x, F0.7, 5x, F0.7)')  self%x, self%y, self%z
    end if

  end subroutine fprintf

  pure subroutine PropZero(self)
    !! 将和与平方和归零
    class(Prop), intent(inout) :: self    !< 输入的Prop类变量

    self%sum = 0.d0
    self%sum2 = 0.d0

  end subroutine PropZero

  pure subroutine PropAccum(self)
    !! 将物理量累加
    class(Prop), intent(inout) :: self    !< 输入的Prop类变量

    self%sum = self%sum + self%val
    self%sum2 = self%sum2 + self%val * self%val

  end subroutine PropAccum

  pure subroutine PropAvg(self, n)
    !! 计算物理量的平均值和平方差
    class(Prop), intent(inout) :: self    !< 输入的Prop类变量
    integer, intent(in) :: n              !< 之前累加的步数

    self%sum = self%sum / n
    self%sum2 = sqrt(Max((self%sum2 / n - self%sum * self%sum), 0.D0))

  end subroutine PropAvg

  pure function SumVectors(Vectors) result(res)
    !! 计算传入的向量数组的各分量和, 并以标量的形式返回
    type(VecR), intent(in), dimension(:) :: Vectors     !< 传入的VecR类数组
    type(VecR) :: res                                   !< 返回的标量VecR类变量

    res%x = sum(Vectors%x)
    res%y = sum(Vectors%y)
    res%z = sum(Vectors%z)

  end function SumVectors

  pure function SumVector(Vector) result(res)
    !! 计算传入的标量的向量各分量之和
    type(VecR), intent(in) :: Vector    !< 传入的VecR类数组
    real(wp) :: res                     !< 返回的各分量之和

    res = Vector%x + Vector%y + Vector%z

  end function SumVector

  pure function SumVec2(Vectors) result(res)
    !! 计算传入的向量数组的各分量的平方和, 并以标量的形式返回
    type(VecR), intent(in), dimension(:) :: Vectors   !< 传入的VecR类数组
    real(wp) :: res                                   !< 返回的向量数组的各分量的平方和

    res = sum(Vectors%x * Vectors%x + Vectors%y * Vectors%y + Vectors%z * Vectors%z)

  end function SumVec2

  pure function VIProd(self) result(res)
    !! 计算VecI的各分量之积
    class(VecI), intent(in) :: self
    integer :: res

    res = self%x * self%y * self%z
  end function VIProd

end module vector
