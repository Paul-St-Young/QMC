double precision function mean(a, n)
  double precision, intent(in) :: a(n)
  integer i
  mean = 0
  do i=1, n
    mean = mean + a(i)
  enddo
  mean = mean/n
end function mean

double precision function stddev(a, n)
  double precision, intent(in) :: a(n)
  double precision mu, mean
  integer i
  mu = mean(a, n)
  stddev = 0
  do i=1, n
    stddev = stddev + (a(i)-mu)**2
  enddo
  stddev = sqrt(stddev/(n-1))
end function stddev
