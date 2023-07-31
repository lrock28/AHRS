# Attitude and Heading Reference System (AHRS)

   * [License](LICENSE)
   * [Changelog](CHANGELOG)

# Description

# Installation

## Matlab
Simply clone or download and the library into your  MATLAB  path folder. MATLAB 2023a was used. Include the dataset if desired for sample data.

# Constants
The following constants are defined from WGS84:

| Description | Variable | Value |
| --- | --- | --- |
| Semi-major axis, m | wgs84.a | 6378137.0 |
| Flattening | wgs84.flattening | 1.0 / 298.257223563 |
| Angular velocity of the Earth, rad/s | wgs84.earth_rot_rate | 7292115.0e-11 |
| Earth's gravitational constant, m^3/s^2 | wgs84.g_const | 3986004.418e8 |

The following constants were also used:

| Description | Variable | Value |
| --- | --- | --- |
| Time Step, seconds | dt_s | 0.01 |
| Magnetic Declination, deg | magvar_deg | -3.26581 |

# Filters

## 6 State EKF AHRS
IMU data was parsed into Matlab from an excel file. After parsing the data, the AHRS attitude was initialized using the accelerometer and magnetometer data. The pitch and roll were calculated from the accelerometer specific force vector, $f_{b}$.

The body acceleration, $a_{b}$, had no way to be calculated, and was thus assumed to be zero. The magnitude of the gravity acceleration vector, $g_n$, was calculated using the current latitude, longitude, and WGS84 ellipsoid data. The yaw was calculated using the three axis magnetometer (TAM) data, $\overrightarrow{h_B}\ $, which needed to be projected onto the horizontal plan using a tilt compensation equation (TC).
        
The roll, $\phi$, and pitch, $\theta$, in the previous equation were used from the accelerometer data.  After the projection was completed, the yaw, $\psi$ was calculated with the projection, $h_{TC}$, and the magnetic declination, $\eta$.        
        
The initial gyroscope bias, $b_{g}$, was given as $b_{g} = [0.00024 -0.00013 .00056]^T rad/s$ based on previous tests with the sensor. The magnetic variation was pulled from the National Oceanic and Atmospheric Association (NOAA) website for the latitude and longitude of the vehicle, given from a GPS receiver. The magnetic variation was assumed to be constant at $-3.2658 rad$, but could vary with the latitude and longitude with a magnetic variation grid.
        
After initializing the INS state, the state was predicted using a time update. The state was comprised of the Euler angles, $\widehat{\psi}$, and the three estimated gyroscope biases, $\widehat{b_g}$. In order to calculate the Euler angles at each time step, the measured gyroscope angular rates, $\overrightarrow{\omega_{B,N\leftarrow}}$, were adjusted for the gyroscope biases and integrated. In order to integrate the angular rates, they were converted to quaternions, discretized at 0.01s intervals, and linearized using a first order Taylor Series. After predicting the Euler angles, the  covariance matrix, $\widehat{P_k}$, could be propagated forward in the time step.
        
The process noise covariance, $Q_k$, is calculated from the power spectral density, S, of the gyroscope, gyroscope characteristics, $\tau_g$ and $\sigma_{bg}$, and the Euler Direction Cosine Matrix, $C_{NB}$.

        
After propagating the covariance matrix forward in time, the AHRS calculates the Euler angles from the accelerometer and magnetometer sensors so as to build the error state and integrate the solutions together.

The innovation, $\delta y_k$, is calculated by taking the residual between the gyroscope Euler angles and the measurement update Euler angles. Then, a Kalman gain, $K$ is calculated in order to correct the original state for the measurement updates using the gyroscope estimated noise covariance, $R$.

With the Kalman gain calculated, a corrected covariance matrix, $P_k$, is then calculated.

Lastly within the measurement update, the attitude and gyroscope bias were updated by multiplying the Kalman gain by the innovation, and adding the error state to the original Euler angles.

# Transformations

## Attitude
These functions convert between:
   * Euler angles ('ZYX' order rotation angles, i.e. yaw, pitch, roll)
   * Direction cosine matrix
   * Quaternions (i, j, k, w)
   * Magnetometer data
   * Accelerometer data

**bfs_euler2quaternion(const T rot1, const T rot2, const T rot3)** Converts a given set of Euler angles to quaternion. Input is the Euler angle vector (order 'ZYX'). Angle input is in radians.

```MATLAB
yaw = 0.7854; 
pitch = 0.1; 
roll = 0;
eul = [yaw, pitch, roll];
quaternion = bfs_euler2quaternion (eul);
```

**bfs_quat2eul(quaternion)** Converts a quaternion to Euler angles. Input is the quaternion. Angles are output in radians.

```MATLAB
q.x = 0.215509;
q.y = 0.432574;
q.z = 0.0846792;
q.w = 0.871358;
quaternion = [q.x, q.y, q.z, q.w];
euler = bfs_quaternion2euler(quaternion);
```

**bfs_quat2dcm(quaternion)** Convert quaternion to direction cosine matrix. Input is the quaternion and output is the direction cosine matrix.

```MATLAB
q.x = 0.215509;
q.y = 0.432574;
q.z = 0.0846792;
q.w = 0.871358;
quaternion = [q.x, q.y, q.z, q.w];
dcm = bfs_quaternion2dcm(quaternion);
```

**tilt(imu)** Convert accelerometer measurements into roll and pitch. Input is accelerometer reading in meters per second. Output is in radians.

```MATLAB
accelmps = [0.1, 0.1, 0.1];
[roll, pitch] = tilt(accelmps);
```

**psi_mag = tilt_comp (hb, phi, theta, magvar_rad)** Convert magnetometer measurements, roll, pitch, and magnetic declination into yaw. Input is magnetometer reading in microtesla and radians. Output is in radians in relation to True North.

```MATLAB
magut = [0.1, 0.1, 0.1];
roll = 0.1;
pitch = 0.1;
magvar_rad = -3.6;
psi = tilt_comp (magut, roll, pitch, magvar_rad);
```
