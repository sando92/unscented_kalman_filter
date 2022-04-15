#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    debug_ = false;

    // false until state x and covariance matrix P are initialized with first measurement received
    is_initialized_ = false;

    // to do or not consistency check
    nis_check_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.5;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 1.5;

    /**
    * DO NOT MODIFY measurement noise values below.
    * These are provided by the sensor manufacturer.
    */

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    /**
    * End DO NOT MODIFY section for measurement noise values 
    */

    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;

    Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ +1);

    // Initialisation of noise covariance matrix
    Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
    Q_(0,0) = std_a_ * std_a_;
    Q_(1,1) = std_yawdd_ * std_yawdd_;

    Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

    weights_ = VectorXd(2 * n_aug_ + 1);
    // set weights
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    double weight = 1 / ( 2 * (lambda_ + n_aug_) );
    weights_(0) = weight_0;
    for (int i=1; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 weights
    weights_(i) = weight;
    }

    n_radar_ = 3;
    // mean predicted measurement
    z_radar_pred_ = VectorXd::Zero(n_radar_);

    R_radar_ = MatrixXd(n_radar_, n_radar_);
    R_radar_ <<  std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;
    // measurement covariance matrix S
    S_radar_ = MatrixXd::Identity(n_radar_, n_radar_);

    // create matrix for sigma points in radar measurement space
    Zsig_radar_ = MatrixXd::Zero(n_radar_, 2 * n_aug_ + 1);

    n_lidar_ = 2;
    // mean predicted measurement
    z_lidar_pred_ = VectorXd::Zero(n_lidar_);

    R_lidar_ = MatrixXd(n_lidar_, n_lidar_);
    R_lidar_ <<  std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;
    // measurement covariance matrix S
    S_lidar_ = MatrixXd::Identity(n_lidar_, n_lidar_);

    // create matrix for sigma points in lidar measurement space
    Zsig_lidar_ = MatrixXd::Zero(n_lidar_, 2 * n_aug_ + 1);

    std::cout << "UKF init done." << std::endl;
}

UKF::~UKF() {}

/**
 * Update UKF from lidar measurement
 */
void UKF::UpdateLidar(const MeasurementPackage &meas_package) {
    if (debug_)
        std::cout << "Update from LIDAR measurement." << std::endl;

    if (is_initialized_) {
        if (use_laser_) {
            // create matrix for cross correlation Tc
            MatrixXd Tc = MatrixXd::Zero(n_x_, n_lidar_);
            Tc.fill(0.0);
            // calculate cross correlation matrix
            for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 sigma points
                // residual
                VectorXd z_diff = Zsig_lidar_.col(i) - z_lidar_pred_;
                NormalizeAngle(&z_diff(1));

                // state difference
                VectorXd x_diff = Xsig_pred_.col(i) - x_;
                NormalizeAngle(&x_diff(3));

                Tc += weights_(i) * x_diff * z_diff.transpose();
            }

            // residual
            VectorXd z_diff = meas_package.raw_measurements_ - z_lidar_pred_;
            NormalizeAngle(&z_diff(1));

            // NIS calculation
            if (nis_check_) {
                double nis = z_diff.transpose() * S_lidar_.inverse() * z_diff;
                if (debug_)
                    std::cout << "NIS is: " << nis << std::endl;
                nis_lidar_.push_back(nis); 
            }

            // Kalman gain K;
            MatrixXd K = Tc * S_lidar_.inverse();

            // update state mean and covariance matrix
            x_ += K * z_diff;
            P_ -= K * S_lidar_ * K.transpose();
        }
    } else {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
        x_(2) = 0.0; // v
        x_(3) = 0.0; // yaw
        x_(4) = 0.0; // yaw_dot

        P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
                0, std_laspy_ * std_laspy_, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;

        std::cout << "State is initialized: " << x_ << std::endl;
    }
}

/**
 * Update UKF from radar measurement.
 */
void UKF::UpdateRadar(const MeasurementPackage &meas_package) {
    if (debug_)
        std::cout << "Update from RADAR measurement." << std::endl;
    if (is_initialized_) {
        if (use_radar_) {
            // create matrix for cross correlation Tc
            MatrixXd Tc = MatrixXd::Zero(n_x_, n_radar_);
            Tc.fill(0.0);
            // calculate cross correlation matrix
            for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 sigma points
                // residual
                VectorXd z_diff = Zsig_radar_.col(i) - z_radar_pred_;
                NormalizeAngle(&z_diff(1));

                // state difference
                VectorXd x_diff = Xsig_pred_.col(i) - x_;
                NormalizeAngle(&x_diff(3));

                Tc += weights_(i) * x_diff * z_diff.transpose();
            }

            // residual
            VectorXd z_diff = meas_package.raw_measurements_ - z_radar_pred_;
            NormalizeAngle(&z_diff(1));

            // NIS calculation
            if (nis_check_) {
                double nis = z_diff.transpose() * S_radar_.inverse() * z_diff;
                if (debug_)
                    std::cout << "NIS is: " << nis << std::endl;
                nis_radar_.push_back(nis);
            }

            // Kalman gain K;
            MatrixXd K = Tc * S_radar_.inverse();

            // update state mean and covariance matrix
            x_ += K * z_diff;
            P_ -= K * S_radar_ * K.transpose();
        }
    } else {
        double rho = meas_package.raw_measurements_(0);
        double phi = meas_package.raw_measurements_(1);
        double rho_dot = meas_package.raw_measurements_(2);

        double p_x = cos( phi ) * rho;
        double p_y = sin( phi ) * rho;

        x_(0) = p_x;
        x_(1) = p_y;
        x_(2) = 0.0; // v
        x_(3) = 0.0; // yaw
        x_(4) = 0.0; // yaw_dot

        P_ << std_radr_* std_radr_, 0, 0, 0, 0,
                0, std_radr_ * std_radr_, 0, 0, 0,
                0, 0, std_radrd_ * std_radrd_, 0, 0,
                0, 0, 0, std_radphi_, 0,
                0, 0, 0, 0, std_radphi_;

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;

        std::cout << "State is initialized: " << x_ << std::endl;
    }
}

/**
 * ProcessMeasurement function
 * update UKF from measurement (RADAR or LIDAR).
 * Once UKF initialized proceed to predict step. 
*/
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
    if ( meas_package.sensor_type_ == MeasurementPackage::LASER )
    {
        if (debug_)
            std::cout << "ProcessMeasurement func called with LASER measurement." << std::endl;
        UpdateLidar(meas_package);
    }
    else if ( meas_package.sensor_type_ == MeasurementPackage::RADAR )
    {
        if (debug_)
            std::cout << "ProcessMeasurement func called with RADAR measurement." << std::endl;
        UpdateRadar(meas_package);
    }

    if (is_initialized_) {
        double dt = (meas_package.timestamp_ - time_us_) * 1e-6;
        time_us_ = meas_package.timestamp_;
        Prediction(dt);
    }
}

/**
 * Create Augmented Sigma points
 */
void UKF::CreateAugmentedSigmaPoints(){
    // create augmented mean state
    VectorXd x_aug = VectorXd::Zero(n_aug_);
    x_aug.head(n_x_) = x_;
    // create augmented covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;

    // create square root matrix
    MatrixXd A = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug_.col(0) = x_aug;
    MatrixXd temp_mat = sqrt(lambda_ + n_aug_) * A;    
    Xsig_aug_.middleCols(1, n_aug_) = temp_mat.colwise() + x_aug;
    Xsig_aug_.middleCols(n_aug_ + 1, n_aug_) = (-temp_mat).colwise() + x_aug;

    if (debug_) {
        std::cout << "CreateAugmentedSigmaPoints done." << std::endl;
        std::cout << "Augmented sigma points are: " << std::endl << Xsig_aug_ << std::endl << std::endl;
    }
}

/**
 * Predict sigma points
 */
void UKF::SigmaPointPrediction(double delta_t){
    // predict sigma points
    for (int i = 0; i < ( 2 * n_aug_ + 1 ); ++i) {
        // extract values for better readability
        double p_x = Xsig_aug_(0,i);
        double p_y = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * ( sin(yaw + yawd * delta_t) - sin(yaw) );
            py_p = p_y + v / yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        // add noise
        px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p += nu_a * delta_t;

        yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p += nu_yawdd * delta_t;

        // write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;

        if (debug_) {    
            std::cout << "SigmaPointPrediction done." << std::endl;
            std::cout << "Predicted sigma points are: " << std::endl << Xsig_pred_ << std::endl << std::endl;
        }
    }
}

/**
 * Predict mean state and covariance matrix from predicted sigma points
 */
void UKF::PredictMeanAndCovariance() {
    // predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x_ +=  weights_(i) * Xsig_pred_.col(i);
    }

    P_.fill(0.0);
    // predicted state covariance matrix
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(&x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }

    if (debug_) {    
        std::cout << "PredictMeanAndCovariance done." << std::endl;
        std::cout << "Predicted state is: " << std::endl << x_ << std::endl;
        std::cout << "Predicted cov matrix is: " << std::endl << P_ << std::endl << std::endl;
    }
}

/**
 * Predict lidar measurement,
 *      Transform sigma points into lidar measurement space
 *      Caculate mean z_lidar_pred_
 *      Calculate corresponding covariance matrix S_lidar_
 */
void UKF::PredictLidarMeasurement() {
    // transform sigma points into measurement space
    Zsig_lidar_.fill(0.0);
    z_lidar_pred_.fill(0.0);
    for (int i = 0; i < (2 * n_aug_ + 1); ++i) {  // 2n+1 sigma points
        // extract values for better readability
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);

        // measurement model
        Zsig_lidar_(0,i) = p_x;
        Zsig_lidar_(1,i) = p_y;
        
        z_lidar_pred_ += weights_(i) * Zsig_lidar_.col(i);
    }

    S_lidar_.fill(0.0);
    for (int i = 0; i < (2 * n_aug_ + 1); ++i) {  // 2n+1 sigma points
        // residual
        VectorXd z_diff = Zsig_lidar_.col(i) - z_lidar_pred_;

        NormalizeAngle(&z_diff(1));

        S_lidar_ += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    S_lidar_ += R_lidar_;
}

/**
 * Predict radar measurement,
 *      Transform sigma points into radar measurement space
 *      Calculate mean z_radar_pred_
 *      Calculate corresponding covariance matrix S_radar_
 */
void UKF::PredictRadarMeasurement() {
    // transform sigma points into measurement space
    Zsig_radar_.fill(0.0);
    z_radar_pred_.fill(0.0);
    for (int i = 0; i < (2 * n_aug_ + 1); ++i) {  // 2n+1 sigma points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                // r
    Zsig_radar_(1,i) = atan2(p_y,p_x);                         // phi
    Zsig_radar_(2,i) = (p_x*v1 + p_y*v2) / Zsig_radar_(0,i);   // r_dot

    z_radar_pred_ += weights_(i) * Zsig_radar_.col(i);
    }

    S_radar_.fill(0.0);
    for (int i = 0; i < (2 * n_aug_ + 1); ++i) {  // 2n+1 sigma points
    // residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_radar_pred_;

    NormalizeAngle(&z_diff(1));

    S_radar_ += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    S_radar_ += R_radar_;
}

/**
 * Main UKF prediction function calls:
 *      CreateAugmentedSigmaPoints()
 *      SigmaPointPrediction()
 *      PredictMeanAndCovariance()
 *      
 *      PredictRadarMeasurement()
 *      PredictLidarMeasurement()
 * 
 */
void UKF::Prediction(double delta_t) {
    if (is_initialized_) {
        if (debug_)
            std::cout << "Prediction time !" << std::endl;

        CreateAugmentedSigmaPoints();
        SigmaPointPrediction(delta_t);

        PredictMeanAndCovariance();

        if (use_radar_)
            PredictRadarMeasurement();
        if (use_laser_)
            PredictLidarMeasurement();
  }
}

/**
 *  Normalize angle between [-PI; PI]
 */
void UKF::NormalizeAngle(double* angle) {
    while (*angle >  M_PI) *angle -= 2. * M_PI;
    while (*angle < -M_PI) *angle += 2. * M_PI;
}