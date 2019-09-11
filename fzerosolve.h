#ifndef __FZEROSOLVE_H__
#define __FZEROSOLVE_H__
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>




//quadratic_root. 
template <class T>
struct quadratic_root_2deriv
{ 

  T c1;
  T c2;


  // Functor returning both 1st and 2nd derivatives.
quadratic_root_2deriv(T const& to_find_root_of, T c2, T c1) : a(to_find_root_of), c1(c1), c2(c2)
  { /* Constructor stores value a to find root of, for example: */ }

  std::tuple<T, T, T> operator()(T const& x)
  {// for example a sphere would be this
  //f(t) 1/4*g^2*t^4-vz*g*t^3+(vx^2+vy^2+vz^2-g*z)*t^2+(2*vx*x+2*vy*y+2*vz*z)*t +x^2 +y^2 +z^2 -R^2   
    //         c4        c3           c2                    c1                            a
    // Return both f(x) and f'(x) and f''(x).
    using namespace boost::math;

    T fx = c2*pow<2>(x)+c1*(x) + a;    // Difference (estimate x^3 - value).
    T dx = 2.*c2*(x)+c1;    // 1st derivative
    T d2x = 2.*c2;  // 2nd derivative
    return std::make_tuple(fx, dx, d2x);  // return fx, dx and d2x.
  }
private:
  T a;                                    // to be solved for. 
};



//find quadratic root
template <class T>
T quadratic_root(T x, T c2, T c1,T guess, T min, T max)
{ 
  // quadratic root
  using namespace std;                  // Help ADL of std functions.
  using namespace boost::math::tools;   // for halley_iterate.



  // Stop when 40% of the digits are correct, can be tuned for your needs of course, hell make it possible to declare for all i care. :
  const int digits = static_cast<int>(std::numeric_limits<T>::digits * 0.4); 
  const boost::uintmax_t maxit = 50;
  boost::uintmax_t it = maxit;
  //magic bit up next. 
  T result = halley_iterate(quadratic_root_2deriv<T>(x,c2,c1), guess, min, max, digits, it);
  return result;
}


//used in cubic_root. 
template <class T>
struct cubic_root_2deriv
{ 
 
  T c1;
  T c2;
  T c3;

  // Functor returning both 1st and 2nd derivatives.
cubic_root_2deriv(T const& to_find_root_of,T c3, T c2, T c1) : a(to_find_root_of), c1(c1), c2(c2), c3(c3)
  { /* Constructor stores value a to find root of, for example: */ }

  std::tuple<T, T, T> operator()(T const& x)
  {// for example a sphere would be this
  //f(t) 1/4*g^2*t^4-vz*g*t^3+(vx^2+vy^2+vz^2-g*z)*t^2+(2*vx*x+2*vy*y+2*vz*z)*t +x^2 +y^2 +z^2 -R^2   
    //         c4        c3           c2                    c1                            a
    // Return both f(x) and f'(x) and f''(x).
    using namespace boost::math;

    T fx = c3*pow<3>(x)+c2*pow<2>(x)+c1*(x) + a;    // Difference (estimate x^3 - value).
    T dx = 3*c3*pow<2>(x)+2*c2*(x)+c1;    // 1st derivative
    T d2x = 6*c3*x+2*c2;  // 2nd derivative
    return std::make_tuple(fx, dx, d2x);  // 'return' fx, dx and d2x.
  }
private:
  T a;                                    // to be solved for. 
};



//find cubic roots, 
template <class T>
T cubic_root(T x, T c3 ,T c2, T c1,T guess, T min, T max)
{ 
  // cubic root
  using namespace std;                  // Help ADL of std functions.
  using namespace boost::math::tools;   // for halley_iterate.



  // Stop when 100% of the digits are correct, can be tuned for your needs of course, hell make it possible to declare for all i care. :
  const int digits = static_cast<int>(std::numeric_limits<T>::digits ); 
  const boost::uintmax_t maxit = 50;
  boost::uintmax_t it = maxit;
  //magic bit up next. 
  T result = halley_iterate(cubic_root_2deriv<T>(x,c3,c2,c1), guess, min, max, digits, it);
  //T result = schroder_iterate(cubic_root_2deriv<T>(x,c3,c2,c1), guess, min, max, digits, it);

  return result;
}



//used in quartic_root. 
template <class T>
struct quartic_root_2deriv
{ 

  T c1;
  T c2;
  T c3;
  T c4;

  // Functor returning both 1st and 2nd derivatives.
quartic_root_2deriv(T const& to_find_root_of, T c4, T c3, T c2, T c1) : a(to_find_root_of), c1(c1), c2(c2), c3(c3), c4(c4)
  { /* Constructor stores value a to find root of, for example: */ }

  std::tuple<T, T, T> operator()(T const& x)
  {// for example a sphere would be this
  //f(t) 1/4*g^2*t^4-vz*g*t^3+(vx^2+vy^2+vz^2-g*z)*t^2+(2*vx*x+2*vy*y+2*vz*z)*t +x^2 +y^2 +z^2 -R^2   
    //         c4        c3           c2                    c1                            a
    // Return both f(x) and f'(x) and f''(x).
    using namespace boost::math;

    T fx = c4*pow<4>(x)+c3*pow<3>(x)+c2*pow<2>(x)+c1*(x) + a;    // Difference (estimate x^3 - value).
    T dx = 4*c4*pow<3>(x) + 3*c3*pow<2>(x)+2*c2*(x)+c1;    // 1st derivative
    T d2x = 12*c4*pow<2>(x)+6*c3*x+2*c2;  // 2nd derivative
    return std::make_tuple(fx, dx, d2x);  // 'return' fx, dx and d2x.
  }
private:
  T a;                                    // to be solved for. 
};



//find quartic root
template <class T>
T quartic_root(T x, T c4, T c3 ,T c2, T c1,T guess, T min, T max)
{ 
  // quartic root
  using namespace std;                  // Help ADL of std functions.
  using namespace boost::math::tools;   // for halley_iterate.



  // Stop when 100% of the digits are correct, can be tuned for your needs of course, hell make it possible to declare for all i care. :
  const int digits = static_cast<int>(std::numeric_limits<T>::digits); 
  const boost::uintmax_t maxit = 50;
  boost::uintmax_t it = maxit;
  //magic bit up next. 
  T result = halley_iterate(quartic_root_2deriv<T>(x,c4,c3,c2,c1), guess, min, max, digits, it);
  return result;
}


#endif