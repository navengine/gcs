#ifndef SATSIM_INCLUDE_GPS_CORRELATOR_SIM_HPP
#define SATSIM_INCLUDE_GPS_CORRELATOR_SIM_HPP 

#include <vector>
#include <complex>
#include <functional>
#include <cassert>
#include <random>

#include <Eigen/Dense>

#include "gps_correlator_model.hpp"


namespace gcs {

template<typename T, int S>
T EvaluatePolynomial(const T arg, const Eigen::Vector<T,S>& coeffs)
{
  T pow = T(1.0);
  T result = T(0.0);
  for (int i = 0; i < S; i++) {
    result += pow * coeffs(i);
    pow *= arg;
  }
  return result;
}

template<typename T, int S>
std::function<T(T)> PolynomialFunction(const Eigen::Vector<T,S>& coeffs)
{
  return std::bind(&EvaluatePolynomial<T,S>, std::placeholders::_1, coeffs);
}


// The first index is defined as the prompt tap, while all others are relative to the prompt
template<typename T = double, typename RandomGen = std::default_random_engine>
class CorrelatorSim
{
public:
  struct Tap
  {
    std::function<T(T)> code_error;
    std::function<T(T)> phase_error;
    std::size_t segments {1};
  };
  
  typedef Eigen::Matrix<std::complex<T>,Eigen::Dynamic,Eigen::Dynamic> CMat;
  typedef Eigen::Vector<std::complex<T>,Eigen::Dynamic> CVec;

private:
  std::vector<Tap> taps_;
  CMat sqrt_cov_;
  bool cov_computed_ {false};
  T corr_period_;
  T cno_;
  
  // noise generation
  std::normal_distribution<T> dist_;
  RandomGen gen_;

  void ComputeCovariance()
  {
    CMat cov = CMat::Identity(taps_.size(), taps_.size());
    for (std::size_t i = 0; i < taps_.size(); i++) {
      for (std::size_t j = 0; j < taps_.size(); j++) {
        if (i == j) continue;
        if (j > i)
          cov(i,j) = PiecewiseI( 
            (taps_[i].segments > taps_[j].segments) ? taps_[i].segments : taps_[j].segments,
            taps_[i].code_error, taps_[j].code_error, taps_[i].phase_error, taps_[j].phase_error,
            corr_period_, cno_
          ) / corr_period_;
        else 
          cov(i,j) = std::conj(cov(j,i));
      }
    }
    std::cout << cov << '\n';
    sqrt_cov_ = cov.llt().matrixL();
    cov_computed_ = true;
  }

public:
  CorrelatorSim() {}

  CorrelatorSim(const std::size_t size) 
  {
    taps_.reserve(size);
  }

  CorrelatorSim(const T corr_period, const T cno)
  {
    corr_period_ = corr_period;
    cno_ = cno;
  }
  
  void SetPeriod(const T corr_period)
  {
    corr_period_ = corr_period;
    cov_computed_ = false;
  }

  // NOTE: C/N0 is in units of Hz, not dB-Hz
  void SetCno(const T cno)
  { cno_ = cno; }

  T Period() const
  { return corr_period_; }

  T Cno() const
  { return cno_; }
  
  void AddTap(const Tap tap)
  {
    taps_.push_back(tap);
    cov_computed_ = false;
  }

  void AddTap(const std::function<T(T)>& code_error,
                     const std::function<T(T)>& phase_error)
  {
    taps_.emplace_back(code_error, phase_error);
    cov_computed_ = false;
  }

  void AddTap(const std::function<T(T)>& code_error,
                     const std::function<T(T)>& phase_error,
                     const std::size_t segments)
  {
    taps_.emplace_back(code_error, phase_error, segments);
    cov_computed_ = false;
  }

  template<int S1, int S2>
  void AddTap(const Eigen::Vector<T,S1>& code_coeffs,
              const Eigen::Vector<T,S2>& carrier_coeffs)
  {
    taps_.emplace_back(PolynomialFunction<T,S1>(code_coeffs),
                       PolynomialFunction<T,S1>(carrier_coeffs));
    cov_computed_ = false;
  }

  template<int S1, int S2>
  void AddTap(const Eigen::Vector<T,S1>& code_coeffs,
              const Eigen::Vector<T,S2>& carrier_coeffs,
              const std::size_t segments)
  {
    taps_.emplace_back(PolynomialFunction<T,S1>(code_coeffs),
                       PolynomialFunction<T,S1>(carrier_coeffs), segments);
    cov_computed_ = false;
  }

  void EraseTap(const std::size_t index)
  {
    assert(index < taps_.size());
    taps_.erase(taps_.begin() + index);
    cov_computed_ = false;
  }

  Tap& GetTap(const std::size_t index)
  { return taps_[index]; }

  void UpdateCovariance()
  { cov_computed_ = false; }

  CVec Simulate()
  {
    if (!cov_computed_) ComputeCovariance();
    CVec noise(taps_.size());
    for (int i = 0; i < taps_.size(); i++) {
      noise(i) = dist_(gen_);
      noise(i) += dist_(gen_) * COMPLEX_I<T>;
    }
    CVec correlations(taps_.size());
    for (int i = 0; i < taps_.size(); i++) {
      correlations(i) = PiecewiseCorrelatorModel(taps_[i].segments, taps_[i].code_error,
                                           taps_[i].phase_error, corr_period_, cno_);
    }
    std::cout << sqrt_cov_ << '\n';
    return correlations + (sqrt_cov_ * noise);
  }

};

} // namespace gcs
#endif
