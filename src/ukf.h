#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage &meas_package);

  /**
   * CreateAugmentedSigmaPoints Creates augmented sigma points
   */
  void CreateAugmentedSigmaPoints();

/**
   * SigmaPointPrediction Predicts sigma points
   * @param delta_t Time between k and k+1 in s
   */
  void SigmaPointPrediction(double delta_t);

    /**
   * PredictMeanAndCovariance Predicts the state, and the state covariance
   * matrix
   */
  void PredictMeanAndCovariance();

  /**
   * PredictRadarMeasurement Predicts radar measurement, tranforms sigma points into
   * radar measurement space to complete the update in a second step.
   */
  void PredictRadarMeasurement();

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage &meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage &meas_package);

    /**
    * Constrains given angle between -pi and pi.
    * @param angle The angle to be constrained
    */
  static double ConstrainAngle(double angle);

  bool debug_;

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // augmented sigma points matrix
  Eigen::MatrixXd Xsig_aug_;

  // noise covariance matrix
  Eigen::MatrixXd Q_;


  // Radar dof
  int n_radar_;
  // mean predicted radar measurement
  Eigen::VectorXd z_radar_pred_;

  // radar noise covariance matrix
  Eigen::MatrixXd R_radar_;
  // radar measurement covariance matrix , innovation covariance matrix Sradar_
  Eigen::MatrixXd S_radar_;

  // sigma points matrix in radar measurement space
  Eigen::MatrixXd Zsig_radar_;


  // Lidar dof
  int n_lidar_;
  // measurement function
  Eigen::MatrixXd H_lidar_;
  // lidar noise covariance matrix
  Eigen::MatrixXd R_lidar_;

  // boolean to choose to do or not consistency check
  bool nis_check_;

  // Containers for consistency check
  std::vector<double> nis_radar_;
  std::vector<double> nis_lidar_;

};

#endif  // UKF_H