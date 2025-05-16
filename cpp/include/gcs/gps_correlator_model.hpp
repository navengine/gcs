#ifndef SATSIM_INCLUDE_GPS_CORR_MODEL_HPP
#define SATSIM_INCLUDE_GPS_CORR_MODEL_HPP

#include <complex>
#include <cmath>
#include <functional>

#include <Eigen/Dense>

namespace gcs {

template <typename T = double>
inline constexpr std::complex<T> COMPLEX_I{
    static_cast<T>(0), static_cast<T>(1)};


// NOTE: this is unit-amplitude
template<typename T>
std::complex<T> StandardGpsCorrModelIntegral(
    const T corr_period,
    const T chip_err,
    const T avg_phase_err,
    const T omega_err)
{
  if (std::abs(chip_err) >= 1.0) return 0.0;
  std::complex<T> result = 1 - std::abs(chip_err);
  if (omega_err != 0.0) {
    T in_sinc = omega_err * corr_period / 2.0;
    result *= std::sin(in_sinc) / in_sinc;
  }
  return result * std::exp(COMPLEX_I<T> * avg_phase_err);
}


template<typename T>
std::complex<T> StandardGpsCorrModel(
      const T corr_period,
      const T chip_err,
      const T avg_phase_err,
      const T omega_err,
      const T cno)
{
  return 2.0 * std::sqrt(cno * corr_period)
    * StandardGpsCorrModelIntegral(corr_period, chip_err, avg_phase_err, omega_err);
}


template<typename T = double>
std::complex<T> CalculateJ(
      const T corr_period,
      const T init_chip_err,
      const T code_rate_err,
      const T init_phase_err,
      const T omega_err)
{
  std::complex<T> phase_exp = std::exp(COMPLEX_I<T> * init_phase_err);
  if (std::abs(omega_err) < 0.0002) {
    if ((init_chip_err + (0.5 * corr_period * code_rate_err)) >= 0.0) {
      return phase_exp * corr_period
        * ( 1.0 - init_chip_err - (0.5 * code_rate_err * corr_period) );
    }
    else {
      return phase_exp * corr_period
        * ( 1.0 + init_chip_err + (0.5 * code_rate_err * corr_period) );
    }
  }
  else {
    if ((init_chip_err + (0.5 * corr_period * code_rate_err)) >= 0.0) {
      std::complex<T> s = COMPLEX_I<T> * omega_err;
      std::complex<T> result = std::exp(s * corr_period)
        * (s * (init_chip_err + (code_rate_err * corr_period) - 1.0) - code_rate_err);
      result -= (s * (init_chip_err - 1.0) - code_rate_err);
      return result * std::exp(COMPLEX_I<T> * init_phase_err)
        / (omega_err * omega_err);
    }
    else {
      std::complex<T> s = COMPLEX_I<T> * omega_err;
      std::complex<T> result = std::exp(s * corr_period)
        * (s * (-init_chip_err - (code_rate_err * corr_period) - 1.0) + code_rate_err);
      result -= (s * (-init_chip_err - 1.0) + code_rate_err);
      return result * std::exp(COMPLEX_I<T> * init_phase_err)
        / (omega_err * omega_err);
    }
  }
}


template<typename T = double>
std::complex<T> CalculateI(
      const T corr_period,
      const T init_chip_err,
      const T code_rate_err,
      const T init_phase_err,
      const T omega_err)
{
  if (code_rate_err == 0.0) {
    return corr_period * StandardGpsCorrModelIntegral(corr_period, init_chip_err,
      init_phase_err + (omega_err * corr_period / 2.0), omega_err);
  }

  T low_lim = (-1.0 - init_chip_err) / code_rate_err;
  T up_lim = (1.0 - init_chip_err) / code_rate_err;
  if (low_lim > up_lim)
    std::swap(low_lim,up_lim);
  if ((low_lim > corr_period) || (up_lim < 0.0)) return 0.0;

  low_lim = (low_lim > 0.0) ? low_lim : 0.0;
  up_lim = (up_lim < corr_period) ? up_lim : corr_period;

  T t0 = -init_chip_err / code_rate_err;

  T chip_offset = init_chip_err + (code_rate_err * low_lim);
  T phase_offset = init_phase_err + (omega_err * low_lim);
  if ( (t0 > low_lim) && (t0 < up_lim) ) {
    std::complex<T> result = CalculateJ(t0 - low_lim, chip_offset, code_rate_err, phase_offset, omega_err);
    chip_offset = init_chip_err + (code_rate_err * t0);
    phase_offset = init_phase_err + (omega_err * t0);
    return result + CalculateJ(up_lim - t0, chip_offset, code_rate_err, phase_offset, omega_err);
  }
  else {
    return CalculateJ(up_lim - low_lim, chip_offset, code_rate_err, phase_offset, omega_err);
  }
}


template<typename T>
std::complex<T> LinearCorrelatorModel(
      const T corr_period,
      const T init_chip_err,
      const T code_rate_err,
      const T init_phase_err,
      const T omega_err,
      const T cno)
{
  return 2.0 * std::sqrt(cno / corr_period)
    * CalculateI(corr_period, init_chip_err, code_rate_err, init_phase_err, omega_err);
}


// code and phase error functions should have input and output of type T
template<typename T = double>
std::complex<T> PiecewiseI(
      const std::size_t segments,
      std::function<T(T)>& code_error, // time [s] -> chips
      std::function<T(T)>& phase_error,  // time [s] -> radians
      const T corr_period,
      const T cno)
{
  std::complex<T> integral = 0.0;
  Eigen::VectorXd times = Eigen::VectorXd::LinSpaced(segments+1,0.0,corr_period);
  for (std::size_t i = 0; i < segments; i++) {
    T dt = times(i+1) - times(i);
    integral += CalculateI(
                dt,
                code_error(times(i)),
                ( code_error(times(i+1)) - code_error(times(i)) ) / dt,
                phase_error(times(i)),
                ( phase_error(times(i+1)) - phase_error(times(i)) ) / dt
      );
  }
  return integral;
}


// version intended for determining the correlation between two signal replicas
template<typename T = double>
std::complex<T> PiecewiseI(
      const std::size_t segments,
      std::function<T(T)>& code1, // time [s] -> chips
      std::function<T(T)>& code2, // time [s] -> chips
      std::function<T(T)>& carrier1,  // time [s] -> radians
      std::function<T(T)>& carrier2,  // time [s] -> radians
      const T corr_period,
      const T cno)
{
  std::complex<T> integral = 0.0;
  Eigen::VectorXd times = Eigen::VectorXd::LinSpaced(segments+1,0.0,corr_period);
  for (std::size_t i = 0; i < segments; i++) {
    T dt = times(i+1) - times(i);
    integral += CalculateI(
                dt,
                code1(times(i)) - code2(times(i)),
                ( code1(times(i+1)) - code2(times(i+1)) - code1(times(i))
                  + code2(times(i)) ) / dt,
                carrier1(times(i)) - carrier2(times(i)),
                ( carrier1(times(i+1)) - carrier2(times(i+1)) 
                  - carrier1(times(i)) + carrier2(times(i)) ) / dt
      );
  }
  return integral;
}


template<typename T = double>
std::complex<T> PiecewiseCorrelatorModel(
      const std::size_t segments,
      std::function<T(T)>& code_error, // time [s] -> chips
      std::function<T(T)>& phase_error,  // time [s] -> radians
      const T corr_period,
      const T cno)
{
  return 2.0 * std::sqrt(cno / corr_period)
         * PiecewiseI(segments,code_error,phase_error,corr_period,cno);
}

/*
template<std::size_t CodeOrder, std::size_t CarrierOrder, typename T = double>
std::complex<T> PiecewiseCorrelatorModel(
                    const std::size_t segments,
                    Polynomial<CodeOrder,1,T>& code_error, // time [s] -> chips
                    Polynomial<CodeOrder,1,T>& phase_error,  // time [s] -> radians
                    const T corr_period,
                    const T cno)
{
  std::complex<T> integral = 0.0;
  Eigen::VectorXd times = Eigen::VectorXd::LinSpaced(segments+1,0.0,corr_period);
  for (std::size_t i = 0; i < segments; i++) {
    T dt = times(i+1) - times(i);
    integral += CalculateI(
                dt,
                code_error(times(i)),
                ( code_error(times(i+1)) - code_error(times(i)) ) / dt,
                phase_error(times(i)),
                ( phase_error(times(i+1)) - phase_error(times(i)) ) / dt
      );
  }
  return 2.0 * std::sqrt(cno / corr_period) * integral;
}
*/

} // namespace gcs
#endif
