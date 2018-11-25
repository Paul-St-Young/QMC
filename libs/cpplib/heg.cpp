class HomogeneousElectronGas
{
public:
  HomogeneousElectronGas()
    : _rs(1.0)
  {}

  double get_rs() const { return _rs; }
  void set_rs(double rs)
  {
    _rs = rs;
  }

private:
  double _rs;
};

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(heg)
{
    class_<HomogeneousElectronGas>("HomogeneousElectronGas")
        .add_property("rs",
          &HomogeneousElectronGas::get_rs,
          &HomogeneousElectronGas::set_rs
        )
    ;
}
